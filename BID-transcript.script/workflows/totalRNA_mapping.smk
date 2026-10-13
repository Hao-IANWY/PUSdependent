import json
import os
import shlex
from snakemake.exceptions import WorkflowError


# 流程的主工作流：
# 1. 从 config.yaml 校验并展开样本定义；
# 2. 完成 trim / 去重 / UMI 区段截取 / 比对；
# 3. 生成各参考类型的 mpileup 表，统一做 gap 校正；
# 4. 针对 T 位点做预筛后，把 Mock/Treat 成对送入 Fisher 检验。
#
# 这里刻意把配置检查放在文件最前面：
# 一旦 config.yaml 的结构不符合约定，希望 Snakemake 在 DAG 构建阶段就直接报错，
# 而不是等到某条 rule 运行时才因为缺字段或路径拼接错误而失败。
def require_mapping(mapping, name):
    if not isinstance(mapping, dict):
        raise WorkflowError(f"Config entry '{name}' must be a mapping.")
    return mapping


def require_string(mapping, key, location):
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise WorkflowError(f"Config entry '{location}.{key}' must be a non-empty string.")
    return value


CONFIG_PATHS = require_mapping(config.get("paths"), "paths")
CONFIG_ENV = require_mapping(config.get("env"), "env")
CONFIG_PIPELINE = require_mapping(config.get("pipeline"), "pipeline")
CONFIG_SAMPLES = require_mapping(config.get("samples"), "samples")

INDEX_BOWTIE2 = require_string(CONFIG_PATHS, "index_bowtie2", "paths")
SCRIPT_DIR = require_string(CONFIG_PATHS, "script_dir", "paths")
PSEUDO_CALLSITES_R = require_string(CONFIG_PATHS, "pseudo_callsites_r", "paths")

PYTHON = require_string(CONFIG_ENV, "python", "env")
R_BIN = require_string(CONFIG_ENV, "r", "env")

T_SITE_WINDOWS = CONFIG_PIPELINE.get("t_site_windows")
WORKDIR = os.getcwd()

FILTER_CONFIG = require_mapping(CONFIG_PIPELINE.get("filter", {}), "pipeline.filter")
MIN_DEPTH = FILTER_CONFIG.get("min_depth", 10)
MIN_AVG_DEPTH = FILTER_CONFIG.get("min_avg_depth", 10)
MAX_PADJ = FILTER_CONFIG.get("max_padj", 0.01)
MIN_GAP = FILTER_CONFIG.get("min_gap", 5)
MIN_TREAT_RATIO = FILTER_CONFIG.get("min_treat_ratio", 0.05)
RATIO_FOLD = FILTER_CONFIG.get("ratio_fold", 2)
MIN_PASS_REPS = FILTER_CONFIG.get("min_pass_reps", 2)

REF_CONFIG = require_mapping(CONFIG_PIPELINE.get("ref"), "pipeline.ref")
REF = {
    ref_type: require_string(REF_CONFIG, ref_type, "pipeline.ref")
    for ref_type in REF_CONFIG
}
REF_TYPES = tuple(REF)

MPILEUP_DIR = "STAT/mpileup/{ref_type}"
ADJUST_GAP_DIR = "STAT/adjustGap/{ref_type}"
FISHER_DIR = "STAT/fisher/{ref_type}"

if not isinstance(T_SITE_WINDOWS, int):
    raise WorkflowError("Config entry 'pipeline.t_site_windows' must be an integer.")

def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


for _key, _value, _label, _positive in (
    ("min_depth", MIN_DEPTH, "positive number", True),
    ("min_avg_depth", MIN_AVG_DEPTH, "positive number", True),
    ("max_padj", MAX_PADJ, "number in (0, 1]", False),
    ("min_gap", MIN_GAP, "non-negative number", False),
    ("min_treat_ratio", MIN_TREAT_RATIO, "non-negative number", False),
    ("ratio_fold", RATIO_FOLD, "positive number", True),
):
    if not _is_number(_value):
        raise WorkflowError(f"Config entry 'pipeline.filter.{_key}' must be a {_label}.")
    if _positive and _value <= 0:
        raise WorkflowError(f"Config entry 'pipeline.filter.{_key}' must be a positive number.")
    if _key == "max_padj" and not (0 < _value <= 1):
        raise WorkflowError("Config entry 'pipeline.filter.max_padj' must be in (0, 1].")
    if _key == "min_gap" and _value < 0:
        raise WorkflowError("Config entry 'pipeline.filter.min_gap' must be a non-negative number.")
    if _key == "min_treat_ratio" and _value < 0:
        raise WorkflowError("Config entry 'pipeline.filter.min_treat_ratio' must be a non-negative number.")

if not isinstance(MIN_PASS_REPS, int) or isinstance(MIN_PASS_REPS, bool) or MIN_PASS_REPS <= 0:
    raise WorkflowError("Config entry 'pipeline.filter.min_pass_reps' must be a positive integer.")

if not REF_TYPES:
    raise WorkflowError("Config entry 'pipeline.ref' must define at least one reference type.")

if len(REF_TYPES) > 2:
    raise WorkflowError(
        "Config entry 'pipeline.ref' must define one or two reference types in mapping order: Bowtie2 first, STAR second."
    )

BOWTIE2_REF_TYPE = REF_TYPES[0]
HAS_STAR_STAGE = len(REF_TYPES) == 2

if HAS_STAR_STAGE:
    STAR_REF_TYPE = REF_TYPES[1]
    INDEX_STAR = require_string(CONFIG_PATHS, "index_star", "paths")
else:
    STAR_REF_TYPE = None
    INDEX_STAR = None


# Preprocessing controls are global for one workflow run. Existing configs remain
# compatible because paired-end input, rmdup, and the historical UMI trimming
# values are used when these keys are omitted.
def parse_read_layout(value):
    if not isinstance(value, str):
        raise WorkflowError("Config entry 'pipeline.read_layout' must be one of: paired, PE, single, SE.")

    aliases = {
        "paired": "paired",
        "pe": "paired",
        "single": "single",
        "se": "single",
    }
    normalized = aliases.get(value.strip().lower())
    if normalized is None:
        raise WorkflowError("Config entry 'pipeline.read_layout' must be one of: paired, PE, single, SE.")
    return normalized


def require_boolean(value, location):
    if not isinstance(value, bool):
        raise WorkflowError(f"Config entry '{location}' must be true or false.")
    return value


def require_nonnegative_integer(value, location):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise WorkflowError(f"Config entry '{location}' must be a non-negative integer.")
    return value


def parse_umi_trim_config(pipeline_config, defaults):
    raw_config = pipeline_config.get("umi_extract")
    if raw_config is None:
        raw_config = {}
    else:
        raw_config = require_mapping(raw_config, "pipeline.umi_extract")

    unsupported_reads = sorted(set(raw_config) - {"R1", "R2"})
    if unsupported_reads:
        raise WorkflowError(
            "Config entry 'pipeline.umi_extract' only supports R1 and R2; "
            f"found: {', '.join(unsupported_reads)}."
        )

    parsed = {}
    for read_name in ("R1", "R2"):
        read_location = f"pipeline.umi_extract.{read_name}"
        read_config = raw_config.get(read_name, {})
        if read_config is None:
            read_config = {}
        else:
            read_config = require_mapping(read_config, read_location)

        unsupported_keys = sorted(set(read_config) - {"trim_front", "trim_tail"})
        if unsupported_keys:
            raise WorkflowError(
                f"Config entry '{read_location}' only supports trim_front and trim_tail; "
                f"found: {', '.join(unsupported_keys)}."
            )

        parsed[read_name] = {
            "trim_front": require_nonnegative_integer(
                read_config.get("trim_front", defaults[read_name]["trim_front"]),
                f"{read_location}.trim_front",
            ),
            "trim_tail": require_nonnegative_integer(
                read_config.get("trim_tail", defaults[read_name]["trim_tail"]),
                f"{read_location}.trim_tail",
            ),
        }

    return parsed


READ_LAYOUT = parse_read_layout(CONFIG_PIPELINE.get("read_layout", "paired"))
RMDUP_ENABLED = require_boolean(CONFIG_PIPELINE.get("rmdup", True), "pipeline.rmdup")
ANALYSIS_READ = "R1" if READ_LAYOUT == "single" else "R2"
DEFAULT_UMI_TRIM = {'R1': {'trim_front': 17, 'trim_tail': 0}, 'R2': {'trim_front': 17, 'trim_tail': 0}}
UMI_TRIM = parse_umi_trim_config(CONFIG_PIPELINE, DEFAULT_UMI_TRIM)
TRIM_GALORE_LENGTH = require_nonnegative_integer(CONFIG_PIPELINE.get("trim_galore_length", 35), "pipeline.trim_galore_length")


# 每个 replicate 的阶段归属（移植自 polyA_mapping.smk）。一个 replicate 可以进入
# 表达量/count 分支、Fisher 位点鉴定分支，或者两者都进；"stages" 键可选，
# 省略时等价于两者都进，因此旧 config 行为完全不变。
#
# 只做 count 的重复（stages: [count]）照常 trim/比对/mpileup/计数/出日志，
# 但绝不会被交给 combine_mpileup_tsv 之后的位点鉴定链 —— 否则加入一个重复会改变
# 其它所有文库的跨样本 gap 校正，迫使整条 Fisher 链全量重跑。
STAGE_COUNT = "count"
STAGE_FISHER = "fisher"
VALID_STAGES = (STAGE_COUNT, STAGE_FISHER)


def parse_replicate_stages(raw_stages, location):
    if raw_stages is None:
        return frozenset(VALID_STAGES)
    if not isinstance(raw_stages, (list, tuple)) or not raw_stages:
        raise WorkflowError(
            f"Config entry '{location}.stages' must be a non-empty list containing "
            + " and/or ".join(VALID_STAGES)
            + "."
        )
    stages = set()
    for stage in raw_stages:
        if not isinstance(stage, str):
            raise WorkflowError(f"Config entry '{location}.stages' entries must be strings.")
        stages.add(stage.strip().lower())
    unsupported = sorted(stages - set(VALID_STAGES))
    if unsupported:
        raise WorkflowError(
            f"Config entry '{location}.stages' only supports {', '.join(VALID_STAGES)}; "
            f"found: {', '.join(unsupported)}."
        )
    return frozenset(stages)


# 把层级化的样本定义摊平成两个索引：
# 1. SAMPLES: 以 file_index 为键，便于单个 replicate 规则直接取输入 fastq；
# 2. FISHER_PAIRS: 以 sample_key 为键，只装进入 Fisher 阶段的重复，保留 Mock/Treat 顺序。
#
# 这里保留 config.yaml 中 replicate 的原始顺序，因为 Pseudo_callsites.R 默认要求：
# Mock 组在前、Treat 组在后，且两个组内部重复次序一一对应。
def build_sample_registry(samples_config):
    sample_records = {}
    fisher_pairs = {}
    count_only_samples = []
    errors = []

    for sample_key, sample_config in samples_config.items():
        if not isinstance(sample_key, str) or not sample_key:
            errors.append("samples keys must be non-empty strings.")
            continue

        sample_location = f"samples.{sample_key}"
        sample_config = require_mapping(sample_config, sample_location)
        sample_id = require_string(sample_config, "sample_id", sample_location)
        if sample_id != sample_key:
            errors.append(
                f"{sample_location}.sample_id must match its key '{sample_key}', got '{sample_id}'."
            )
            continue

        fisher_pair = {}
        all_groups = {}
        replicate_names = {}

        for treatment_key, treatment_name in (("mock", "Mock"), ("treat", "Treat")):
            treatment_location = f"{sample_location}.{treatment_key}"
            replicate_config = require_mapping(sample_config.get(treatment_key), treatment_location)
            if not replicate_config:
                errors.append(f"{treatment_location} must define at least one replicate.")
                continue

            group_members = []
            fisher_members = []
            treatment_replicates = []

            for replicate_name, replicate_info in replicate_config.items():
                replicate_location = f"{treatment_location}.{replicate_name}"
                replicate_info = require_mapping(replicate_info, replicate_location)
                file_index = require_string(replicate_info, "file_index", replicate_location)
                read1 = require_string(replicate_info, "read1", replicate_location)
                stages = parse_replicate_stages(replicate_info.get("stages"), replicate_location)
                if READ_LAYOUT == "paired":
                    read2 = require_string(replicate_info, "read2", replicate_location)
                else:
                    read2 = replicate_info.get("read2")
                    if read2 is not None and (not isinstance(read2, str) or not read2):
                        raise WorkflowError(
                            f"Config entry '{replicate_location}.read2' must be omitted or be a non-empty string."
                        )

                if file_index in sample_records:
                    errors.append(f"Duplicate file_index '{file_index}' found in {replicate_location}.")
                    continue

                sample_records[file_index] = {
                    "sample_key": sample_key,
                    "treatment": treatment_name,
                    "replicate": replicate_name,
                    "read1": read1,
                    "read2": read2,
                    "stages": stages,
                }
                group_members.append(file_index)
                if STAGE_FISHER in stages:
                    fisher_members.append(file_index)
                treatment_replicates.append(replicate_name)

            all_groups[treatment_name] = tuple(group_members)
            fisher_pair[treatment_name] = tuple(fisher_members)
            replicate_names[treatment_name] = tuple(treatment_replicates)

        mock_files = all_groups.get("Mock", ())
        treat_files = all_groups.get("Treat", ())
        mock_replicates = replicate_names.get("Mock", ())
        treat_replicates = replicate_names.get("Treat", ())

        if not mock_files or not treat_files:
            errors.append(f"{sample_location} must define both mock and treat replicates.")
            continue
        if len(mock_files) != len(treat_files):
            errors.append(
                f"{sample_location} has {len(mock_files)} mock replicates but {len(treat_files)} treat replicates."
            )
            continue
        if mock_replicates != treat_replicates:
            errors.append(
                f"{sample_location} replicate names must match between mock and treat; got {list(mock_replicates)} vs {list(treat_replicates)}."
            )
            continue

        mock_fisher = fisher_pair.get("Mock", ())
        treat_fisher = fisher_pair.get("Treat", ())

        if not mock_fisher and not treat_fisher:
            # 该样本的所有重复都是 count-only：整个样本退出 Fisher 的 DAG。
            count_only_samples.append(sample_key)
            continue
        if not mock_fisher or not treat_fisher:
            missing_side = "mock" if not mock_fisher else "treat"
            errors.append(
                f"{sample_location}: '{missing_side}' has no replicate marked for the {STAGE_FISHER} "
                f"stage. Mark at least one replicate on each side, or mark the whole sample "
                f"count-only with 'stages: [{STAGE_COUNT}]'."
            )
            continue

        fisher_pairs[sample_key] = fisher_pair

    if errors:
        raise WorkflowError("Invalid config.yaml sample definition:\n- " + "\n- ".join(errors))

    if count_only_samples:
        print(
            f"[totalRNA] sample(s) excluded from the {STAGE_FISHER} branch (all replicates are "
            f"{STAGE_COUNT}-only): " + ", ".join(count_only_samples)
        )

    return sample_records, fisher_pairs


SAMPLES, FISHER_PAIRS = build_sample_registry(CONFIG_SAMPLES)
FILE_INDEX = tuple(sorted(SAMPLES))
FISHER_SAMPLE_KEYS = tuple(FISHER_PAIRS)
FISHER_FILE_INDEX = tuple(
    sorted(
        {
            file_index
            for pair in FISHER_PAIRS.values()
            for group_members in pair.values()
            for file_index in group_members
        }
    )
)
COUNT_ONLY_FILE_INDEX = tuple(
    file_index for file_index in FILE_INDEX if file_index not in set(FISHER_FILE_INDEX)
)

# 位点鉴定链（combine_mpileup_tsv -> adjust_gap -> filter_t_sites -> fisher_sites）
# 一律消费 FISHER_FILE_INDEX，绝不碰 FILE_INDEX：因此新增一个 count-only 文库只会
# 多出比对/计数/日志 job，COMBINED 表、adjustGap 的样本集合与全部 Fisher 产物保持不变。


# Optional raw expression counting. The feature is intentionally off unless a
# project sets pipeline.expression.enabled, so legacy configs keep their exact
# previous rule-all targets. Two matrices are emitted per reference type: the
# Mock matrix keeps the historical mock-only sample layout, while all_samples
# covers every configured library. Sample order follows config.yaml insertion
# order through SAMPLES above.
EXPRESSION_GROUP_LABELS = ("Mock", "all_samples")


def parse_expression_config(pipeline_config, samples, ref_types):
    expression_config = pipeline_config.get("expression", {})
    if expression_config is None:
        expression_config = {}
    else:
        expression_config = require_mapping(
            expression_config, "pipeline.expression"
        )

    enabled = expression_config.get("enabled", False)
    if not require_boolean(enabled, "pipeline.expression.enabled"):
        return False, {}, {}

    mock_file_indexes = tuple(
        file_index
        for file_index, sample_info in samples.items()
        if sample_info["treatment"].lower() == "mock"
    )
    if not mock_file_indexes:
        raise WorkflowError(
            "pipeline.expression requires at least one mock sample for the "
            "Mock count matrix; check the samples section."
        )
    group_file_indexes = {
        "Mock": mock_file_indexes,
        "all_samples": tuple(samples),
    }

    count_mode_config = require_mapping(
        expression_config.get("count_mode"), "pipeline.expression.count_mode"
    )
    unsupported_keys = tuple(
        sorted((set(count_mode_config) - set(ref_types)), key=str)
    )
    missing_keys = tuple(sorted((set(ref_types) - set(count_mode_config)), key=str))
    if unsupported_keys or missing_keys:
        details = []
        if unsupported_keys:
            details.append(f"unsupported: {', '.join(map(str, unsupported_keys))}")
        if missing_keys:
            details.append(f"missing: {', '.join(map(str, missing_keys))}")
        raise WorkflowError(
            "pipeline.expression.count_mode keys must exactly match pipeline.ref; "
            + "; ".join(details)
        )

    count_modes = {
        ref_type: require_string(
            count_mode_config, ref_type, "pipeline.expression.count_mode"
        )
        for ref_type in ref_types
    }
    invalid_modes = tuple(
        sorted(
            (ref_type for ref_type, mode in count_modes.items() if mode not in VALID_EXPRESSION_MODES),
            key=str,
        )
    )
    if invalid_modes:
        raise WorkflowError(
            "pipeline.expression.count_mode values must be one of: "
            f"{', '.join(VALID_EXPRESSION_MODES)}. Invalid reference types: "
            f"{', '.join(invalid_modes)}"
        )

    return True, count_modes, group_file_indexes


VALID_EXPRESSION_MODES = ("all_alignments", "primary")
(
    EXPRESSION_ENABLED,
    EXPRESSION_COUNT_MODES,
    EXPRESSION_GROUP_FILE_INDEXES,
) = parse_expression_config(CONFIG_PIPELINE, SAMPLES, REF_TYPES)
EXPRESSION_COUNT_OUTPUTS = tuple(
    f"count/{ref_type}/{group}.{EXPRESSION_COUNT_MODES[ref_type]}_counts.tsv"
    for ref_type in REF_TYPES
    for group in EXPRESSION_GROUP_LABELS
) if EXPRESSION_ENABLED else ()


def expression_count_output(ref_type, group):
    if ref_type not in EXPRESSION_COUNT_MODES:
        raise WorkflowError(
            f"No expression count mode configured for reference type '{ref_type}'."
        )
    if group not in EXPRESSION_GROUP_FILE_INDEXES:
        raise WorkflowError(f"Unknown expression count group '{group}'.")
    return (
        f"count/{ref_type}/{group}."
        f"{EXPRESSION_COUNT_MODES[ref_type]}_counts.tsv"
    )


def validate_expression_count_request(wildcards):
    expected_output = expression_count_output(wildcards.ref_type, wildcards.group)
    observed_output = (
        f"count/{wildcards.ref_type}/{wildcards.group}."
        f"{wildcards.mode}_counts.tsv"
    )
    if observed_output != expected_output:
        raise WorkflowError(
            f"Requested expression output '{observed_output}' does not match "
            f"the configured output '{expected_output}'."
        )


def expression_count_bams(wildcards):
    validate_expression_count_request(wildcards)
    return [
        f"BAM/{file_index}.{wildcards.ref_type}.realign.filtered.BAM"
        for file_index in EXPRESSION_GROUP_FILE_INDEXES[wildcards.group]
    ]


def expression_count_bam_args(wildcards):
    bam_inputs = expression_count_bams(wildcards)
    file_indexes = EXPRESSION_GROUP_FILE_INDEXES[wildcards.group]
    return " ".join(
        f"--bam {shlex.quote(f'{file_index}={bam_path}')}"
        for file_index, bam_path in zip(file_indexes, bam_inputs)
    )


def normalize_raw_path(path_value):
    if path_value.startswith("/"):
        return path_value
    return os.path.abspath(os.path.join(WORKDIR, path_value))



def raw_read_path(wildcards, read_name):
    read_key = {"R1": "read1", "R2": "read2"}[read_name]
    read_path = SAMPLES[wildcards.file_index].get(read_key)
    if not read_path:
        raise WorkflowError(
            f"Sample '{wildcards.file_index}' does not define {read_key}, "
            f"which is required for {READ_LAYOUT}-end processing."
        )
    return normalize_raw_path(read_path)


def clean_fastq_path(file_index, read_name):
    return f"fix.fastq/{file_index}.{read_name}.clean.fq.gz"


def dedup_fastq_path(file_index, read_name):
    return f"fix.fastq/{file_index}.{read_name}.dup.fq.gz"


def pre_umi_fastq_path(file_index, read_name):
    if RMDUP_ENABLED:
        return dedup_fastq_path(file_index, read_name)
    return clean_fastq_path(file_index, read_name)


def umi_fastq_path(file_index, read_name):
    return f"umi_extract.fastq/{file_index}.{read_name}.extract.fq.gz"

# 通用路径模板。后面的辅助函数都围绕这些模板展开，避免各 rule 自己重复拼路径。
FIX_FASTQ_PREFIX = "fix.fastq/{file_index}"
UMI_FASTQ_PREFIX = "umi_extract.fastq/{file_index}"
BAM_PREFIX = "BAM/{file_index}"
BAM_REF_PREFIX = "BAM/{file_index}.{ref_type}"
MPILEUP_PREFIX = "STAT/mpileup/{ref_type}/{file_index}.{ref_type}"
LOG_PREFIX = "logs/{file_index}"
DELETION_SUMMARY = "logs/deletion_summary.tsv"
PROCESSING_SUMMARY = "logs/processing_summary.tsv"


def ref_path(wildcards):
    return REF[wildcards.ref_type]


def render_path(template, wildcards):
    return template.format(**wildcards)


def build_mpileup_tsv_path(file_index, ref_type, outdir=None):
    base_dir = outdir if outdir is not None else MPILEUP_DIR.format(ref_type=ref_type)
    return f"{base_dir}/{file_index}.{ref_type}.mpileup.tsv"


def build_filter_t_path(file_index, ref_type, outdir=None):
    base_dir = outdir if outdir is not None else ADJUST_GAP_DIR.format(ref_type=ref_type)
    return f"{base_dir}/{file_index}.{ref_type}.filter_T.tsv"


# 同一个 ref_type 下的所有 replicate 会先汇总为一个 COMBINED 表，
# 后续 adjustGap 用它统一估计位点层面的 gap 校正结果。
def get_mpileup_tsv_inputs(wildcards):
    mpileup_dir = render_path(MPILEUP_DIR, wildcards)
    return [build_mpileup_tsv_path(file_index, wildcards.ref_type, mpileup_dir) for file_index in FISHER_FILE_INDEX]


# 以下一组函数为 logs / processing_stats / deletion_stats 三条规则准备输入：
# 把每个文库（以及每个 ref_type）的日志路径打包成脚本需要的 JSON 参数，
# 避免在 rule 里重复展开 FILE_INDEX x REF_TYPES。
def deletion_log_dataset_records():
    records = []
    for file_index in FILE_INDEX:
        sample = SAMPLES[file_index]
        for ref_type in REF_TYPES:
            records.append(
                {
                    "sample": sample["treatment"],
                    "replicate": sample["replicate"],
                    "file_index": file_index,
                    "ref_type": ref_type,
                    "log": LOG_PREFIX.format(file_index=file_index) + ".logs.tsv",
                }
            )
    return records


def processing_log_dataset_records():
    records = []
    for file_index in FILE_INDEX:
        sample = SAMPLES[file_index]
        records.append(
            {
                "sample": sample["treatment"],
                "replicate": sample["replicate"],
                "file_index": file_index,
                "log": LOG_PREFIX.format(file_index=file_index) + ".logs.tsv",
            }
        )
    return records


def get_all_log_inputs(wildcards):
    return [
        LOG_PREFIX.format(file_index=file_index) + ".logs.tsv"
        for file_index in FILE_INDEX
    ]


def build_deletion_log_dataset_args(wildcards):
    return " ".join(
        "--log-dataset " + shlex.quote(json.dumps(record, sort_keys=True, separators=(",", ":")))
        for record in deletion_log_dataset_records()
    )


def build_processing_log_dataset_args(wildcards):
    return " ".join(
        "--dataset " + shlex.quote(json.dumps(record, sort_keys=True, separators=(",", ":")))
        for record in processing_log_dataset_records()
    )


def build_library_bam_args(wildcards):
    args = []
    for ref_type in REF_TYPES:
        for stage in ("forward", "realign", "filtered"):
            value = f"{ref_type},{stage},BAM/{wildcards.file_index}.{ref_type}.{stage}.bam"
            if stage == "filtered":
                value = f"{ref_type},{stage},BAM/{wildcards.file_index}.{ref_type}.realign.filtered.BAM"
            args.append("--bam " + shlex.quote(value))
    return " ".join(args)


def build_library_mpileup_args(wildcards):
    return " ".join(
        "--mpileup " + shlex.quote(f"{ref_type},{build_mpileup_tsv_path(wildcards.file_index, ref_type)}")
        for ref_type in REF_TYPES
    )


def build_fisher_cli_args(sample_key, ref_type):
    adjust_gap_dir = ADJUST_GAP_DIR.format(ref_type=ref_type)
    sample = FISHER_PAIRS[sample_key]
    input_groups = []
    name_groups = []

    for treatment in ("Mock", "Treat"):
        members = sample[treatment]
        input_groups.append(":".join(build_filter_t_path(file_index, ref_type, adjust_gap_dir) for file_index in members))
        name_groups.append(":".join(members))

    return ",".join(input_groups), ",".join(name_groups)


# Pseudo_callsites.R 要求输入按照 “Mock 组全部重复 -> Treat 组全部重复” 的顺序传入。
# 这里直接复用 build_fisher_cli_args 的展开逻辑，避免 rule.input 和 rule.params 各写一套。
def get_fisher_inputs(wildcards):
    input_files, _ = build_fisher_cli_args(wildcards.sample_key, wildcards.ref_type)
    return [path for group in input_files.split(",") for path in group.split(":")]


def get_fisher_input_files(wildcards):
    input_files, _ = build_fisher_cli_args(wildcards.sample_key, wildcards.ref_type)
    return input_files


def get_fisher_input_names(wildcards):
    _, input_names = build_fisher_cli_args(wildcards.sample_key, wildcards.ref_type)
    return input_names


# 最终目标覆盖两部分：
# 1. 每个参考类型的中间统计表（便于排查流程某一步是否异常）；
# 2. 每个 sample_key 的 Fisher 结果和 mapping 日志统计。
rule all:
    input:
        expand(
            f"{MPILEUP_PREFIX}.mpileup.tsv",
            file_index=FILE_INDEX,
            ref_type=REF_TYPES,
        ),
        expand(
            f"{MPILEUP_DIR}/COMBINED-{{ref_type}}.tsv",
            ref_type=REF_TYPES,
        ),
        expand(
            f"{ADJUST_GAP_DIR}/.done",
            ref_type=REF_TYPES,
        ),
        expand(
            f"{ADJUST_GAP_DIR}/{{file_index}}.{{ref_type}}.filter_T.tsv",
            file_index=FISHER_FILE_INDEX,
            ref_type=REF_TYPES,
        ),
        expand(
            f"{FISHER_DIR}/{{sample_key}}_sites.Fisher.tsv",
            ref_type=REF_TYPES,
            sample_key=FISHER_SAMPLE_KEYS,
        ),
        expand(
            f"{FISHER_DIR}/{{sample_key}}_sites.positive.tsv",
            ref_type=REF_TYPES,
            sample_key=FISHER_SAMPLE_KEYS,
        ),
        expand(
            f"{LOG_PREFIX}.logs.tsv",
            file_index=FILE_INDEX,
        ),
        PROCESSING_SUMMARY,
        DELETION_SUMMARY,
        *EXPRESSION_COUNT_OUTPUTS
# totalRNA/tRNA 文库的有效信息在 R2：8bp UMI + TATAGGG + insert + polyA + adapter。
# 读取、修剪与比对阶段只关心单个 replicate 的处理，因此全部以 file_index 为主键。
# totalRNA/tRNA libraries are analyzed from R2 for PE data. In SE mode the
# sole read is provided as read1 and follows the same trimming/mapping path.
if READ_LAYOUT == "paired":
    rule trim:
        input:
            raw_fq=lambda wildcards: raw_read_path(wildcards, "R2"),
        output:
            clean_fq=f"{FIX_FASTQ_PREFIX}.R2.clean.fq.gz",
        threads: 8
        shell:
            """
            trim_galore --cores {threads} --length {TRIM_GALORE_LENGTH} --adapter AAAAAAAAAAA \
                --basename {wildcards.file_index} \
                -o fix.fastq/ \
                --fastqc {input.raw_fq}
            mv fix.fastq/{wildcards.file_index}_trimmed.fq.gz {output.clean_fq}
            """
else:
    rule trim:
        input:
            raw_fq=lambda wildcards: raw_read_path(wildcards, "R1"),
        output:
            clean_fq=f"{FIX_FASTQ_PREFIX}.R1.clean.fq.gz",
        threads: 8
        shell:
            """
            trim_galore --cores {threads} --length {TRIM_GALORE_LENGTH} --adapter AAAAAAAAAAA \
                --basename {wildcards.file_index} \
                -o fix.fastq/ \
                --fastqc {input.raw_fq}
            mv fix.fastq/{wildcards.file_index}_trimmed.fq.gz {output.clean_fq}
            """


if RMDUP_ENABLED:
    rule rmdup:
        input:
            fix_fq=lambda wildcards: clean_fastq_path(wildcards.file_index, ANALYSIS_READ),
        output:
            dup_fq=f"{FIX_FASTQ_PREFIX}.{ANALYSIS_READ}.dup.fq.gz",
        threads: 8
        shell:
            """
            seqkit rmdup {input.fix_fq} -s -j {threads} -o {output.dup_fq}
            """


rule umi_extract:
    input:
        fq=lambda wildcards: pre_umi_fastq_path(wildcards.file_index, ANALYSIS_READ),
    output:
        extract_fq=f"{UMI_FASTQ_PREFIX}.{ANALYSIS_READ}.extract.fq.gz",
    params:
        trim_front=lambda wildcards: UMI_TRIM[ANALYSIS_READ]["trim_front"],
        trim_tail=lambda wildcards: UMI_TRIM[ANALYSIS_READ]["trim_tail"],
    threads: 8
    shell:
        """
        fastp --thread {threads} -AQL --html /dev/null --json /dev/null \
            --in1 {input.fq} \
            --out1 {output.extract_fq} \
            --trim_front1 {params.trim_front} \
            --trim_tail1 {params.trim_tail}
        """


rule bowtie2_map:
    input:
        lambda wildcards: umi_fastq_path(wildcards.file_index, ANALYSIS_READ)
    output:
        bam=f"{BAM_PREFIX}.{BOWTIE2_REF_TYPE}.bam",
        un=f"{BAM_PREFIX}.{BOWTIE2_REF_TYPE}.unmap.fq",
        report=f"{BAM_PREFIX}.{BOWTIE2_REF_TYPE}.report",
    threads: 16
    shell:
        """
        export LC_ALL=C
        bowtie2 -p {threads} \
            --local --ma 2 --score-min G,10,7 -D 20 -R 3 -L 8 -N 1 -i S,1,0.5 --mp 6,3 --rdg 1,2 --rfg 6,3 --norc -a \
            --no-unal \
            --un {output.un} \
            -x {INDEX_BOWTIE2} \
            -U {input} 2>{output.report} | \
            {SCRIPT_DIR}/samFilter | \
            samtools view -O BAM -o {output.bam}
        """


if HAS_STAR_STAGE:
    rule mrna_map:
        input:
            f"{BAM_PREFIX}.{BOWTIE2_REF_TYPE}.unmap.fq"
        output:
            bam=f"{BAM_PREFIX}.{STAR_REF_TYPE}.bam",
            log=f"{BAM_PREFIX}.{STAR_REF_TYPE}.Log.final.out",
        threads: 24
        shell:
            """
            STAR --runThreadN {threads} --genomeDir {INDEX_STAR} \
                --alignEndsType EndToEnd \
                --scoreDelOpen -1 \
                --scoreDelBase -1 \
                --scoreInsOpen -2 \
                --scoreInsBase -2 \
                --outFilterMatchNmin 15 \
                --outFilterMatchNminOverLread 0.8 \
                --outFilterMismatchNmax 10 \
                --outFilterMismatchNoverLmax 0.2 \
                --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                --alignSJDBoverhangMin 1 \
                --alignSJoverhangMin 5 \
                --chimSegmentMin 20 \
                --chimOutType WithinBAM HardClip \
                --chimJunctionOverhangMin 15 \
                --chimScoreJunctionNonGTAG 0 \
                --outFilterMultimapNmax 10 \
                --outFilterMultimapScoreRange 0 \
                --outSAMmultNmax -1 \
                --outMultimapperOrder Random \
                --outReadsUnmapped Fastx \
                --outSAMtype BAM Unsorted \
                --outStd BAM_Unsorted \
                --outSAMattrRGline ID:{wildcards.file_index} SM:{wildcards.file_index} LB:RNA PL:Illumina PU:SE \
                --outSAMattributes NH HI AS nM NM MD jM jI MC ch \
                --readFilesIn {input} \
                --outFileNamePrefix BAM/{wildcards.file_index}.{STAR_REF_TYPE}. \
                > {output.bam}
            """


# 由于后续是在转录组参考上寻找正链 T 位点，且当前文库是链特异性的，
# 这里明确丢弃 SAM flag 16 的比对，避免把反向链 reads 带入 gap 统计。
rule filter_reverse:
    input:
        bam=f"{BAM_REF_PREFIX}.bam"
    output:
        f"{BAM_REF_PREFIX}.forward.bam"
    threads: 4
    shell:
        """
        samtools view -@ {threads} -b -F 16 {input.bam} > {output}
        """


rule gap_realign:
    input:
        bam=f"{BAM_REF_PREFIX}.forward.bam"
    output:
        f"{BAM_REF_PREFIX}.realign.bam"
    params:
        ref=ref_path
    shell:
        """
        {SCRIPT_DIR}/realignGap -r {params.ref} -i {input.bam} -o {output}
        """


rule sort_cal_filter_bam:
    input:
        bam=f"{BAM_REF_PREFIX}.realign.bam"
    output:
        f"{BAM_REF_PREFIX}.realign.filtered.BAM"
    params:
        ref=ref_path
    threads: 4
    shell:
        """
        samtools sort -@ {threads} -m 4G {input.bam} | \
            samtools calmd -@ {threads} - {params.ref} 2>/dev/null | \
            samtools view -@ {threads} \
                --reference {params.ref} \
                -e '[NM]<=5 && [NM]/(qlen-sclen)<=0.1' \
                -O BAM \
                -o {output}
        """


rule bam_index:
    input:
        bam=f"{BAM_REF_PREFIX}.realign.filtered.BAM"
    output:
        bai=f"{BAM_REF_PREFIX}.realign.filtered.BAM.bai"
    threads: 4
    shell:
        """
        samtools index -@ {threads} {input.bam} {output.bai}
        """


rule bam_mpileup:
    input:
        bam=f"{BAM_REF_PREFIX}.realign.filtered.BAM",
        bai=f"{BAM_REF_PREFIX}.realign.filtered.BAM.bai"
    output:
        mpileup=f"{MPILEUP_PREFIX}.mpileup",
        tsv=f"{MPILEUP_PREFIX}.mpileup.tsv",
    params:
        ref=ref_path,
        outdir=lambda wildcards: render_path(MPILEUP_DIR, wildcards)
    shell:
        """
        mkdir -p {params.outdir}
        samtools mpileup --reverse-del -d 0 -BQ0 -f {params.ref} {input.bam} -o {output.mpileup} --ff UNMAP,QCFAIL -aa
        awk '$3 ~ /^[ACGTNacgtn]$/ {{print}}' {output.mpileup} \
            | {SCRIPT_DIR}/cpup -H -s -i \
            | awk '$4 == "+"' > {output.tsv}
        """


# 每个参考类型只执行一次 COMBINED 聚合：
# 这样 adjustGap 在同一个 ref_type 下看到的是完整样本集合，而不是单个 replicate 的局部统计。
rule combine_mpileup_tsv:
    input:
        get_mpileup_tsv_inputs
    output:
        combined=f"{MPILEUP_DIR}/COMBINED-{{ref_type}}.tsv"
    params:
        outdir=lambda wildcards: render_path(MPILEUP_DIR, wildcards),
        suffix=lambda wildcards: f".{wildcards.ref_type}.mpileup.tsv",
    shell:
        """
        mkdir -p {params.outdir}
        {PYTHON} {SCRIPT_DIR}/combine_mpileup_tsv.py \
            --output {output.combined} \
            --suffix '{params.suffix}' \
            {input}
        """


# 候选位点收紧链路从这里开始：
# combine_mpileup_tsv -> adjust_gap -> filter_t_sites -> fisher_sites。
# adjustGap 会为每个 replicate 写出独立 adjusted.tsv；`.done` 只是一个额外的完成标记，
# 用于让 `rule all` 能清楚表达“这一层的批量校正已经跑完”。
rule adjust_gap:
    input:
        combined=f"{MPILEUP_DIR}/COMBINED-{{ref_type}}.tsv"
    output:
        adjusted=expand(
            f"{ADJUST_GAP_DIR}/{{file_index}}.adjusted.tsv",
            file_index=FISHER_FILE_INDEX,
            ref_type="{ref_type}",
        ),
        done=f"{ADJUST_GAP_DIR}/.done",
    params:
        outdir=lambda wildcards: render_path(ADJUST_GAP_DIR, wildcards)
    shell:
        """
        mkdir -p {params.outdir}
        {SCRIPT_DIR}/adjustGap -i {input.combined} -o {params.outdir}
        touch {output.done}
        """


# 单个 replicate 的 adjusted.tsv 会在这里进一步规整成 Fisher 输入格式：
# 1. 只保留参考碱基为 T 且 depth >= pipeline.filter.min_depth 的记录；
# 2. 以位点为中心提取上下文序列；
# 3. 产出标准列顺序 chr/pos/strand/context/ratio/depth/gap。
#
# 这一步不做组间判断，只负责把每个 replicate 自身的可检验位点准备好。
rule filter_t_sites:
    input:
        adjusted=f"{ADJUST_GAP_DIR}/{{file_index}}.adjusted.tsv"
    output:
        filtered=f"{ADJUST_GAP_DIR}/{{file_index}}.{{ref_type}}.filter_T.tsv"
    params:
        ref_genome=ref_path,
        windows=T_SITE_WINDOWS,
        min_depth=MIN_DEPTH
    shell:
        """
        set -euo pipefail
        output_dir=$(dirname "{output.filtered}")
        mkdir -p "$output_dir"
        tmp_dir=$(mktemp -d "$output_dir/.filter_t_sites.XXXXXX")
        trap 'rm -rf "$tmp_dir"' EXIT

        # 先只保留目标 T 位点，并要求 depth >= min_depth。
        awk 'NR>1 && $4=="T" && $6>={params.min_depth}' "{input.adjusted}" > "$tmp_dir/t_sites.tsv"

        # 以位点为中心构造上下文窗口，后面用参考序列补上 context 列。
        awk -v w="{params.windows}" 'BEGIN{{OFS="\t"}} {{
            start=$2-w-1; if (start<0) start=0;
            end=$2+w;
            id=$1"_"$2"_"$3;
            print $1,start,end,id,0,$3
        }}' "$tmp_dir/t_sites.tsv" > "$tmp_dir/context_windows.bed"

        bedtools slop -i "$tmp_dir/context_windows.bed" -g "{params.ref_genome}.fai" -b 0 > "$tmp_dir/context_windows.fixed.bed"
        bedtools getfasta -fi "{params.ref_genome}" -bed "$tmp_dir/context_windows.fixed.bed" -s -bedOut > "$tmp_dir/contexts.bed"
        awk 'BEGIN{{OFS="\t"}} {{print $0,$1"_"$2"_"$3}}' "$tmp_dir/t_sites.tsv" > "$tmp_dir/t_sites.with_id.tsv"

        # 按位点 id 把提取到的 context 拼回原表，生成后续 Fisher 使用的标准输入。
        {{
            printf 'chr\tpos\tstrand\tcontext\tratio\tdepth\tgap\n'
            csvtk join -t -f "8;4" "$tmp_dir/t_sites.with_id.tsv" "$tmp_dir/contexts.bed" | \
                awk -F '\t' 'BEGIN{{OFS="\t"}} {{
                    $4=$NF;
                    NF--;
                    if (NF > 7) NF=7;
                    print
                }}'
        }} > "{output.filtered}"
        """


# Fisher 检验以 sample_key 为单位聚合同一生物学样本的所有重复。
# 这里显式保留 Mock -> Treat 的顺序契约，因为 Pseudo_callsites.R 会按这个顺序切组。
# 另外需要注意：`positive.tsv` 不是命令行单独指定的输出，而是 Pseudo_callsites.R
# 根据 `-o *.Fisher.tsv` 自动派生出的副产物，所以 output 中两个文件都属于脚本接口约定的一部分。
rule fisher_sites:
    input:
        get_fisher_inputs
    output:
        fisher=f"{FISHER_DIR}/{{sample_key}}_sites.Fisher.tsv",
        positive=f"{FISHER_DIR}/{{sample_key}}_sites.positive.tsv"
    params:
        outdir=lambda wildcards: render_path(FISHER_DIR, wildcards),
        input_files=get_fisher_input_files,
        input_names=get_fisher_input_names,
        min_avg_depth=MIN_AVG_DEPTH,
        max_padj=MAX_PADJ,
        min_gap=MIN_GAP,
        min_treat_ratio=MIN_TREAT_RATIO,
        ratio_fold=RATIO_FOLD,
        min_pass_reps=MIN_PASS_REPS,
    shell:
        """
        mkdir -p {params.outdir}
        {R_BIN} --vanilla -f {PSEUDO_CALLSITES_R} --args \
            -i "{params.input_files}" \
            -n "{params.input_names}" \
            --min-avg-depth {params.min_avg_depth} \
            --max-padj {params.max_padj} \
            --min-gap {params.min_gap} \
            --min-treat-ratio {params.min_treat_ratio} \
            --ratio-fold {params.ratio_fold} \
            --min-pass-reps {params.min_pass_reps} \
            -o {output.fisher}
        """


# BAM/mpileup 阶段的路径都按 REF_TYPES 顺序展开，
# 这样 logs 规则可以按位置把它们 zip 回各参考类型的名字上。
def stage_bam_paths(file_index, suffix):
    return [f"BAM/{file_index}.{ref_type}.{suffix}" for ref_type in REF_TYPES]


def stage_mpileup_paths(file_index):
    return [build_mpileup_tsv_path(file_index, ref_type) for ref_type in REF_TYPES]


# Two raw expression matrices per reference type: Mock keeps the historical
# mock-only layout, all_samples spans every configured library.
# all_alignments preserves the historical ncRNA/rRNA counting behavior, while
# primary prevents a multimapping read from contributing to multiple mRNAs.
rule expression_counts:
    input:
        fasta=lambda wildcards: REF[wildcards.ref_type],
        bams=expression_count_bams
    output:
        counts="count/{ref_type}/{group}.{mode}_counts.tsv"
    params:
        mode=lambda wildcards: EXPRESSION_COUNT_MODES[wildcards.ref_type],
        bam_args=expression_count_bam_args
    threads: 16
    shell:
        """
        set -euo pipefail
        {PYTHON} {SCRIPT_DIR}/count_feature_alignments.py \
            --fasta {input.fasta} \
            --mode {params.mode} \
            --jobs {threads} \
            --output {output.counts} \
            {params.bam_args}
        """


# 该规则只负责把 mapping 各阶段的读段数与比对统计汇总成一张表，
# 方便后续 QC 或比较不同样本批次；它不参与候选位点筛选本身。
# 除各阶段读数外，build_library_log.py 还会基于 filtered BAM 与 mpileup TSV
# 统计每个 ref_type 的缺失（deletion/gap）信号：
#   <ref>_deletion_signal_reads/_ratio        带缺失的 reads 数及其占比
#   <ref>_assessed_bases/_gap_signal_bases    有覆盖的位点数 / 其中带 gap 信号的位点数
#   <ref>_gap_signal_total/_base_ratio        gap 信号总数 / 位点层面缺失率
rule logs:
    input:
        script=f"{SCRIPT_DIR}/build_library_log.py",
        raw_fq=lambda wildcards: raw_read_path(wildcards, ANALYSIS_READ),
        trimmed_fq=lambda wildcards: clean_fastq_path(wildcards.file_index, ANALYSIS_READ),
        rm_dup_reads=lambda wildcards: pre_umi_fastq_path(wildcards.file_index, ANALYSIS_READ),
        bowtie2_log=f"{BAM_PREFIX}.{BOWTIE2_REF_TYPE}.report",
        star_log=(
            f"{BAM_PREFIX}.{STAR_REF_TYPE}.Log.final.out"
            if HAS_STAR_STAGE
            else []
        ),
        forward_bams=lambda wildcards: stage_bam_paths(wildcards.file_index, "forward.bam"),
        realign_bams=lambda wildcards: stage_bam_paths(wildcards.file_index, "realign.bam"),
        filtered_bams=lambda wildcards: stage_bam_paths(wildcards.file_index, "realign.filtered.BAM"),
        mpileup_tsvs=lambda wildcards: stage_mpileup_paths(wildcards.file_index),
    output:
        f"{LOG_PREFIX}.logs.tsv"
    threads: 4
    params:
        builder=f"{SCRIPT_DIR}/build_library_log.py",
        sample=lambda wildcards: SAMPLES[wildcards.file_index]["treatment"],
        replicate=lambda wildcards: SAMPLES[wildcards.file_index]["replicate"],
        read_layout=READ_LAYOUT,
        rmdup=str(RMDUP_ENABLED).lower(),
        bam_args=build_library_bam_args,
        mpileup_args=build_library_mpileup_args,
        star_arg=(
            lambda wildcards: (
                f"--star-log BAM/{wildcards.file_index}.{STAR_REF_TYPE}.Log.final.out"
                if HAS_STAR_STAGE
                else ""
            )
        )
    shell:
        """
        set -euo pipefail
        {PYTHON} {params.builder} \
            --output {output} \
            --sample {params.sample} \
            --replicate {params.replicate} \
            --file-index {wildcards.file_index} \
            --read-layout {params.read_layout} \
            --rmdup {params.rmdup} \
            --raw-fq {input.raw_fq} \
            --trimmed-fq {input.trimmed_fq} \
            --dedup-fq {input.rm_dup_reads} \
            --bowtie2-log {input.bowtie2_log} \
            {params.star_arg} \
            --samtools-threads {threads} \
            {params.bam_args} \
            {params.mpileup_args}
        """


# 把所有文库的 logs.tsv 按阶段汇总为一张长表（Sample/Replicate/Stage/Reads），
# 用来检查 raw -> trimmed -> dedup -> 各参考比对 各环节 reads 的损耗情况。
rule processing_stats:
    input:
        script=f"{SCRIPT_DIR}/summarize_processing_logs.py",
        logs=get_all_log_inputs
    output:
        PROCESSING_SUMMARY
    params:
        summarizer=f"{SCRIPT_DIR}/summarize_processing_logs.py",
        datasets=build_processing_log_dataset_args,
        ref_types=",".join(REF_TYPES)
    shell:
        """
        set -euo pipefail
        {PYTHON} {params.summarizer} \
            --ref-types {params.ref_types} \
            --output {output} \
            {params.datasets}
        """


# 从各文库 logs.tsv 中抽取缺失信号相关列，按 file_index x ref_type 汇总成一张表，
# 并复核缺失 reads 数 / 缺失位点数不超过各自的分母。
rule deletion_stats:
    input:
        script=f"{SCRIPT_DIR}/summarize_deletion_signals.py",
        logs=get_all_log_inputs
    output:
        DELETION_SUMMARY
    params:
        summarizer=f"{SCRIPT_DIR}/summarize_deletion_signals.py",
        datasets=build_deletion_log_dataset_args
    shell:
        """
        set -euo pipefail
        {PYTHON} {params.summarizer} \
            --output {output} \
            {params.datasets}
        """
