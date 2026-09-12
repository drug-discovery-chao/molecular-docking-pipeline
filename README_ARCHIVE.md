# Archive - 早期版本说明

> 此目录存放项目早期（v1）的独立脚本，保留用于展示从初学者到工程化项目的演进过程。

---

## v1 特征（2026年8月）

### 两个独立脚本
- `mtor_drug_screening_v1.py` — mTOR 靶点
- `sirt1_drug_screening_v1.py` — SIRT1 靶点

### v1 能力
- ChEMBL API 数据拉取
- RDKit 计算 8 种分子描述符（MW/LogP/HBD/HBA/RotatableBonds/TPSA/NumRings/NumAromRings）
- Lipinski 五规则 + Veber 规则双重过滤
- 简单线性扣分评分（0-100）
- 输出 CSV 结果

### v1 不足（v2 改进方向）
1. **代码重复**：两个脚本 90% 相同，仅靶点不同
2. **去重逻辑粗糙**：只保留首次出现的 pChEMBL，丢失更优活性值
3. **无异常处理**：ChEMBL API 超时会导致整个脚本崩溃
4. **无日志系统**：全用 print，生产环境不可用
5. **无 PAINS 过滤**：可能混入假阳性干扰化合物
6. **评分模型简单**：线性扣分，无量纲加权
7. **无可视化**：结果纯表格，无法直观展示
8. **无多靶点对比**：各跑各的，无法横向比较

---

## v2 改进（`../drug_screening_pipeline.py`）

详见根目录主脚本，改进点：
- 统一为 `DrugScreeningPipeline` 类，支持任意靶点
- PAINS 过滤（假阳性去除）
- QED（定量类药性）+ SA Score（合成可及性）
- 6 维度加权评分模型
- 异常处理 + logging
- 5 种可视化图表
- 多靶点对比分析

---

*保留早期版本是为了让评审者看到代码能力和项目思维的演进。*
