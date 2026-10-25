# trua4_dependent_trna_heatmap — trua4 dependent 叶绿体 tRNA 结构热图

## 目的

只使用 leaf 组织中 `trua4` dependent 的叶绿体 tRNA 位点，展示这些位点在 Pt tRNA
二级结构布局上的分布（行 = Pt tRNA 家族，列 = canonical 0–76 位，
值 = WT − Mutant 的修饰比例差）。

对应论文版 `trua4_dependent_trna_heatmap`。

## 输入

- `tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_1-leaf.dependent-sites.Sprinzl.tsv`
- `tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_2-leaf.dependent-sites.Sprinzl.tsv`
  （已由 `add_sprinzl_pos.py` 加上 `Sprinzl_pos`）
- `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`（Pt 家族代表，表达量口径）

## 处理规则

1. 两条株系分别读取 dependent-sites 表，过滤 `is_dependent == TRUE` 且 `chr ^tRNA-Pt-`。
2. 按 `chr + pos + strand` 合并两株系：
   - `Mut_Ratio` = 可用株系 `Mut_Ratio` 的均值；
   - `dependency_support` = `shared` / `trua4_1_only` / `trua4_2_only`；
   - `Delta_Ratio(WT-Mut) = WT_Ratio − Mut_Ratio`。
3. **fix_pos 直接取上游 `Sprinzl_pos`**：不再在本目录做坐标转换，也不依赖结构表；
   两株系按 `chr + pos + strand` 合并时 `fix_pos = coalesce(fix_pos_x, fix_pos_y)`。
4. 位点按 `chr → family` 归到 `family_top` 的 Pt 家族（Pt 为 1 家族 : 1 拷贝）；
   同一行同一 canonical 列有多个位点时取最大 Delta。
5. 色标范围自适应：`ramp_max = max(0.3, ceiling(max(Delta) * 10) / 10)`，
   避免高 Delta 位点全部饱和；当前数据为 0–0.5。

## 当前结果

7 个 dependent 位点（6 shared + 1 `trua4_2_only`）落在 7 个 Pt tRNA 家族的反密码子茎
canonical 38–40，Delta 范围 0.107–0.465；7 个位点的 `Sprinzl_pos` 与旧版结构表映射逐一相同。
行标签改为家族名（如 `Pt-Asn-GTT`，不再带拷贝号 `-1`）。

## 输出

- `trua4_pt_trna_dependent_union.tsv`：两株系合并后的 dependent 位点
  （含 `fix_pos`、每株系 Mut_Ratio、`dependency_support`）
- `trua4_pt_trna_dependent_mapped_sites.tsv`：匹配到 Pt 代表结构后的位点
  （含 `Rep_ID`、`display_label`、`canonical_pos`、`matrix_col_id`）
- `trua4_pt_trna_dependent_heatmap.pdf/png/svg`（18 in 宽，高度随行数自适应）

热图底部注释含区段色带、canonical 位号与区段名；右侧独立带边框的竖直色标。
77 列在 18 in 宽度下较密，区段名用 14 pt、位号用 12 pt，避免相邻区段（尤其
Anticodon Stem / Loop / Stem）压字。

## 运行

```bash
bash data.sh
```
