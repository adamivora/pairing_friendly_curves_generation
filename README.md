# Pairing-friendly Curves Generation

This repository tries to implement (in SageMath) the methods from the standard [ISO/IEC 15946-5:2017 - Elliptic curve generation](https://www.iso.org/standard/69726.html), Chapter 7: _Constructing elliptic curves by complex multiplication_.

The current state of the implementation is following:

|                                                Method | Status |
| ----------------------------------------------------: | :----: |
| [General construction](WIP/complex_multiplication.py) |   🚧   |
|   [Miyaji-Nakabayashi-Takano (MNT) curve](WIP/mnt.py) |   🚧   |
|      [Barreto-Naehrig (BN) curve](barreto_naehrig.py) |   ✅   |
|                               Freeman curve (F curve) |   ❌   |
|              [Cocks-Pinch (CP) curve](cocks_pinch.py) |   ✅   |

✅ = completed, 🚧 = work-in-progress, ❌ = not started

## How to run

For example, you can run the Barreto-Naehrig method as follows:

```bash
sage barreto_naehrig.py -h
```

This command prints the help message which shows all the available parameters of the generation method.
