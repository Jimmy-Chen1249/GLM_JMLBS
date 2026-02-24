## Folder structure

### `Code_Clayton_One`
Simulation and estimation code for the joint model where the dependence between two survival times is modeled via a **Clayton copula**.  
In this setting, the shared random effect \(b_i\) is **1-dimensional**.

### `Code_Clayton_Two`
Extension of `Code_Clayton_One` to the case where the shared random effect \(b_i\) is **2-dimensional** (e.g., random intercept + random slope).  
The dependence between the two survival times is still modeled via a **Clayton copula**.

### `Code_robust` (Robust)
Robustness study under the **1-dimensional** shared random effect \(b_i\).  
We investigate robustness across three copula families for modeling dependence between the two survival times:
- **Gaussian copula**
- **Clayton copula**
- **FGM copula**
