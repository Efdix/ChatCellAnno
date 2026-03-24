from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QFileDialog,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QMessageBox,
    QLineEdit,
    QWidget,
    QScrollArea,
)
import os


class SourceFileRow(QWidget):
    """单个一级文件配置行：别名 + 路径 + 浏览 + 删除。"""

    def __init__(self, parent_plugin, index):
        super().__init__()
        self.parent_plugin = parent_plugin
        self.index = index

        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)

        self.alias_edit = QLineEdit(f"一级文件{index}")
        self.alias_edit.setPlaceholderText("一级文件别名")

        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText("选择一级文件 CSV...")

        self.browse_btn = QPushButton("浏览")
        self.browse_btn.clicked.connect(self.browse_file)

        self.remove_btn = QPushButton("移除")
        self.remove_btn.clicked.connect(self.remove_self)

        row.addWidget(QLabel("名称:"))
        row.addWidget(self.alias_edit)
        row.addWidget(QLabel("路径:"))
        row.addWidget(self.path_edit)
        row.addWidget(self.browse_btn)
        row.addWidget(self.remove_btn)

    def browse_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "选择一级文件", "", "CSV Files (*.csv)")
        if path:
            self.path_edit.setText(path)

    def remove_self(self):
        self.parent_plugin.remove_source_row(self)


class CellOriginTracePlugin:
    """
    通用细胞来源追溯插件。

    目标：
    - 支持 1..N 个一级文件追溯到一个目标文件中的某个 cluster。
    - 去除物种特定命名，改为“一级文件/目标文件”语义。
    - 支持规范格式并兼容历史格式。
    """

    SUPPORTED_CLUSTER_HINTS = (
        "cluster_id",
        "cluster",
        "leiden",
        "seurat_clusters",
        "res",
    )

    SUPPORTED_CELL_ID_HINTS = (
        "cell_id",
        "cell",
        "barcode",
        "barcodes",
    )

    def __init__(self, main_window):
        self.main_window = main_window
        self.source_rows = []

    def init_ui(self):
        if hasattr(self.main_window, "gear_menu"):
            action = self.main_window.gear_menu.addAction("🔍 Cell Origin Tracer (Generic)")
            action.triggered.connect(self.show_dialog)

    def show_dialog(self):
        self.dialog = QDialog(self.main_window)
        self.dialog.setWindowTitle("Cell Origin Tracer (通用追溯)")
        self.dialog.resize(1050, 700)

        root = QVBoxLayout(self.dialog)

        hint = QLabel(
            "建议规范格式: 列名包含 cell_id 与 cluster_id。\n"
            "兼容历史格式: 允许第一列/索引为细胞ID，聚类列可命名为 leiden/cluster 等。"
        )
        hint.setStyleSheet("color: #444; padding: 4px;")
        root.addWidget(hint)

        # 目标文件
        target_row = QHBoxLayout()
        self.target_path_edit = QLineEdit()
        self.target_path_edit.setPlaceholderText("选择目标文件（例如整合结果或子集结果）")
        btn_target = QPushButton("浏览目标文件")
        btn_target.clicked.connect(lambda: self.browse_file_to(self.target_path_edit, "选择目标文件"))
        target_row.addWidget(QLabel("目标文件:"))
        target_row.addWidget(self.target_path_edit)
        target_row.addWidget(btn_target)
        root.addLayout(target_row)

        # 查询 cluster
        cluster_row = QHBoxLayout()
        self.target_cluster_input = QLineEdit()
        self.target_cluster_input.setPlaceholderText("输入目标文件中的 Cluster ID（例如: 0）")
        cluster_row.addWidget(QLabel("目标 Cluster:"))
        cluster_row.addWidget(self.target_cluster_input)
        root.addLayout(cluster_row)

        # 一级文件区域
        src_header = QHBoxLayout()
        src_header.addWidget(QLabel("一级文件列表:"))
        self.btn_add_source = QPushButton("+ 添加一级文件")
        self.btn_add_source.clicked.connect(self.add_source_row)
        src_header.addWidget(self.btn_add_source)
        src_header.addStretch(1)
        root.addLayout(src_header)

        self.sources_container = QWidget()
        self.sources_layout = QVBoxLayout(self.sources_container)
        self.sources_layout.setContentsMargins(0, 0, 0, 0)
        self.sources_layout.setSpacing(6)

        self.sources_scroll = QScrollArea()
        self.sources_scroll.setWidgetResizable(True)
        self.sources_scroll.setWidget(self.sources_container)
        self.sources_scroll.setMinimumHeight(170)
        root.addWidget(self.sources_scroll)

        # 默认一个一级文件
        self.add_source_row()

        # 运行按钮
        run_row = QHBoxLayout()
        self.btn_run = QPushButton("开始追溯")
        self.btn_run.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")
        self.btn_run.clicked.connect(self.run_trace)
        run_row.addWidget(self.btn_run)
        run_row.addStretch(1)
        root.addLayout(run_row)

        # 结果表
        self.result_table = QTableWidget()
        self.result_table.setColumnCount(7)
        self.result_table.setHorizontalHeaderLabels(
            [
                "一级文件",
                "原始Cluster",
                "目标Cluster细胞数",
                "占目标Cluster比例(%)",
                "该原始Cluster总细胞数",
                "目标细胞覆盖该原始Cluster比例(%)",
                "匹配细胞示例",
            ]
        )
        self.result_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        root.addWidget(self.result_table)

        self.dialog.exec()

    def browse_file_to(self, edit_widget, title):
        path, _ = QFileDialog.getOpenFileName(self.main_window, title, "", "CSV Files (*.csv)")
        if path:
            edit_widget.setText(path)

    def add_source_row(self):
        row = SourceFileRow(self, len(self.source_rows) + 1)
        self.source_rows.append(row)
        self.sources_layout.addWidget(row)

    def remove_source_row(self, row_widget):
        if len(self.source_rows) <= 1:
            QMessageBox.warning(self.dialog, "提示", "至少需要保留一个一级文件。")
            return
        self.sources_layout.removeWidget(row_widget)
        self.source_rows = [r for r in self.source_rows if r is not row_widget]
        row_widget.setParent(None)
        row_widget.deleteLater()
        # 重新编号默认名
        for i, row in enumerate(self.source_rows, start=1):
            row.index = i
            if not row.alias_edit.text().strip() or row.alias_edit.text().startswith("一级文件"):
                row.alias_edit.setText(f"一级文件{i}")

    def _normalize_df(self, df, file_label):
        """
        规范化 DataFrame 到统一字段：
        - cell_id 列
        - cluster_col（返回列名）
        """
        # 兼容把 cell_id 放在索引的情况
        if df.index.name is not None and str(df.index.name).strip().lower() in self.SUPPORTED_CELL_ID_HINTS:
            df = df.reset_index()

        # 尝试定位 cell_id 列
        cell_col = None
        for col in df.columns:
            low = str(col).strip().lower()
            if low in self.SUPPORTED_CELL_ID_HINTS:
                cell_col = col
                break

        # 兼容历史格式：第一列通常是 cell barcode
        if cell_col is None and len(df.columns) > 0:
            first_col = df.columns[0]
            if str(first_col).startswith("Unnamed") or str(first_col).strip() == "":
                cell_col = first_col
            else:
                # 若第一列唯一值接近行数，也视作 cell_id
                uniq_ratio = df[first_col].nunique(dropna=True) / max(1, len(df))
                if uniq_ratio > 0.95:
                    cell_col = first_col

        if cell_col is None:
            raise ValueError(f"{file_label}: 未识别到 cell_id 列，请添加 'cell_id' 列。")

        # 尝试定位 cluster 列
        cluster_col = None
        for col in df.columns:
            low = str(col).strip().lower()
            if low in self.SUPPORTED_CLUSTER_HINTS:
                cluster_col = col
                break
        if cluster_col is None:
            for col in df.columns:
                low = str(col).strip().lower()
                if "cluster" in low or "leiden" in low or "res" in low:
                    cluster_col = col
                    break

        if cluster_col is None:
            raise ValueError(f"{file_label}: 未识别到 cluster 列，请添加 'cluster_id' 列。")

        norm = df.copy()
        norm["__cell_id__"] = norm[cell_col].astype(str)
        norm["__cluster__"] = norm[cluster_col].astype(str)

        return norm, cluster_col

    def run_trace(self):
        try:
            import pandas as pd
        except Exception:
            QMessageBox.critical(self.dialog, "缺少依赖", "当前环境缺少 pandas，请先安装：pip install pandas")
            return

        target_path = self.target_path_edit.text().strip()
        target_cluster = self.target_cluster_input.text().strip()

        if not target_path or not os.path.exists(target_path):
            QMessageBox.warning(self.dialog, "错误", "请先选择有效的目标文件。")
            return
        if not target_cluster:
            QMessageBox.warning(self.dialog, "错误", "请输入目标 Cluster。")
            return

        # 收集一级文件配置
        source_specs = []
        for row in self.source_rows:
            alias = row.alias_edit.text().strip() or f"一级文件{row.index}"
            path = row.path_edit.text().strip()
            if not path:
                continue
            source_specs.append((alias, path))

        if not source_specs:
            QMessageBox.warning(self.dialog, "错误", "请至少提供一个一级文件。")
            return

        # 读取并标准化目标文件
        try:
            df_target = pd.read_csv(target_path)
            norm_target, target_cluster_col = self._normalize_df(df_target, "目标文件")
        except Exception as e:
            QMessageBox.critical(self.dialog, "错误", f"读取目标文件失败: {e}")
            return

        target_cells_df = norm_target[norm_target["__cluster__"] == str(target_cluster)]
        target_cells = set(target_cells_df["__cell_id__"].tolist())
        total_target = len(target_cells)

        if total_target == 0:
            QMessageBox.information(
                self.dialog,
                "提示",
                f"目标文件列 '{target_cluster_col}' 中未找到 Cluster={target_cluster} 的细胞。",
            )
            return

        results = []

        for alias, source_path in source_specs:
            if not os.path.exists(source_path):
                results.append([alias, "[文件不存在]", 0, 0.0, 0, 0.0, "-"])
                continue

            try:
                df_src = pd.read_csv(source_path)
                norm_src, src_cluster_col = self._normalize_df(df_src, alias)
            except Exception as e:
                results.append([alias, f"[格式错误] {e}", 0, 0.0, 0, 0.0, "-"])
                continue

            matched = norm_src[norm_src["__cell_id__"].isin(target_cells)]
            if matched.empty:
                results.append([alias, "[无匹配细胞]", 0, 0.0, 0, 0.0, "-"])
                continue

            cluster_counts = matched["__cluster__"].value_counts()
            all_cluster_counts = norm_src["__cluster__"].value_counts()

            for src_cluster, cnt in cluster_counts.items():
                src_total = int(all_cluster_counts.get(src_cluster, 0))
                in_target_pct = round(cnt / total_target * 100, 2)
                source_coverage_pct = round(cnt / src_total * 100, 2) if src_total > 0 else 0.0

                sample_cells = matched[matched["__cluster__"] == src_cluster]["__cell_id__"].head(3).tolist()
                sample_text = ", ".join(sample_cells)

                results.append(
                    [
                        alias,
                        str(src_cluster),
                        int(cnt),
                        in_target_pct,
                        src_total,
                        source_coverage_pct,
                        sample_text if sample_text else "-",
                    ]
                )

        # 排序：按目标贡献比例降序
        results.sort(key=lambda x: x[3], reverse=True)

        self.result_table.setRowCount(len(results))
        for i, row in enumerate(results):
            for j, value in enumerate(row):
                self.result_table.setItem(i, j, QTableWidgetItem(str(value)))

        QMessageBox.information(
            self.dialog,
            "完成",
            f"追溯完成。目标Cluster={target_cluster}，目标细胞数={total_target}。",
        )


def register_plugin(main_window):
    return CellOriginTracePlugin(main_window)
