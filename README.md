# BABY-1L: the run-6 tSIE background, applied to runs 3 and 4

The one-page write-up is [studies/release_identifiability/README.md](studies/release_identifiability/README.md). Fits, controls and caveats are in [APPENDIX.md](studies/release_identifiability/APPENDIX.md).

To reproduce from the public LIBRA data (about 10 minutes on a laptop):

```
cd studies/release_identifiability
make venv && make data && make reproduce-fast
```

`make reproduce-full` reruns everything at the replicate counts quoted in the write-up (about 2 h on 64 cores).

The raw LSC data, scripts and blank scan come from the [LIBRA-project](https://github.com/LIBRA-project) BABY-1L repositories and [libra-toolbox](https://github.com/LIBRA-project/libra-toolbox) (MIT). The paper is [arXiv 2509.26174](https://arxiv.org/abs/2509.26174).

Prepared by Mike German with AI assistance (Claude).
