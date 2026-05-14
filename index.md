# HGNChelper: identification and correction of invalid gene symbols for human and mouse

Please cite our software:

Oh S, Abdelnabi J, Al-Dulaimi R et al. HGNChelper: identification and
correction of invalid gene symbols for human and mouse. F1000Research
2020, 9:1493 (<https://doi.org/10.12688/f1000research.28033.1>)

## Updating

To update the symbols maps for human and mouse yourself, download this
repository and run:

`./update.sh`

from its root directory. Note that this script uses the “roxygen2” R
library to update the documentation.

Alternatively, you can use updated maps without updating the package,
see
[`?getCurrentMaps`](https://waldronlab.io/HGNChelper/reference/getCurrentMaps.md).
