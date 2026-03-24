# Cell Origin Tracer 输入格式规范

本规范用于通用来源追溯插件（Cell Origin Tracer）。

## 1. 推荐规范（强烈建议）
每个 CSV 文件至少包含两列：

- `cell_id`: 细胞唯一标识（字符串）
- `cluster_id`: 聚类标签（字符串或数字）

示例：

```csv
cell_id,cluster_id,annotation
cell_0001,0,Mesenchyme
cell_0002,3,Mesenchyme
```

## 2. 兼容格式（历史数据）
插件也兼容以下情况：

- 细胞 ID 在第一列（例如 `Unnamed: 0`）
- 聚类列命名为：
  - `leiden_res.3.0`
  - `cluster`
  - `leiden`
  - `seurat_clusters`
  - 或其他包含 `cluster` / `leiden` / `res` 关键词的列名

## 3. 追溯逻辑说明

- 目标文件：用于选择目标 Cluster（例如整合结果或子集结果）
- 一级文件：一个或多个来源文件
- 插件会用 `cell_id` 交集匹配目标 Cluster 的细胞在各一级文件中的来源分布

输出包含两类比例：

1. `占目标Cluster比例(%)`
   - 公式：`目标Cluster中来自该原始Cluster的细胞数 / 目标Cluster细胞总数`

2. `目标细胞覆盖该原始Cluster比例(%)`
   - 公式：`目标Cluster中来自该原始Cluster的细胞数 / 该原始Cluster在一级文件中的总细胞数`

## 4. 数据清洗建议

- 确保 `cell_id` 在不同文件中使用同一命名体系（不要混用前缀/后缀）。
- 建议先去重：同一文件中 `cell_id` 不应重复。
- 建议统一 cluster 类型为字符串，以减少类型比较歧义。
