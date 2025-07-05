1. Rewrite Pretrain/setup.py as follows.
```python
import sys
sys.path.append('/content/CAFE-MPP')
```

2. Rewrite Prediction/setup.py as follows.
```python
# Prediction/setup.py
from pathlib import Path
from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np

here = Path(__file__).parent.resolve()

extensions = [
    Extension(
        "Prediction.algos",                    # ★ 完全修飾モジュール名を明示
        [str(here / "algos.pyx")],             #   相対パスではなく絶対パスで安全に
        include_dirs=[np.get_include()],       #   NumPy ヘッダー
    )
]

setup(
    name="Prediction-algos",
    ext_modules=cythonize(
        extensions,
        language_level="3",                    # Py3 構文を確実に
        annotate=False,                        # .html が要るなら True に
    ),
    zip_safe=False,
)

```

3. Rewrite Prediction/trainer.py as follows

```python
import sys
sys.path.append('/content/CAFE-MPP')
```

```
from collections.abc import Iterable
# from collections import Iterable
```

