# ctsCausal

### Overview

ctsCausal is a tool for nonparametric estimation and inference of a continuous dose-response curve, as well as for null hypothesis testing of whether a dose-response is flat or not. It has two main functions ``ctsCausal`` and ``ctsCausalTest``.

The function for nonparametric estimation and inference, ``ctsCausal``, is implemented with three different methods. The method ``dr.isoreg`` is for [causal isotonic regression](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12372), the method ``dr.loclin`` is for [doubly robust estimation of continuous treatment effects](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12212) and the method ``dr.debiased`` is for debiased doubly robust estimation of continuous treatment effects.

``ctsCausalTest`` has two different methods for null hypothesis testing. ``dr.nd`` is for [nonparametric tests of causal null with nondiscrete exposures](https://www.tandfonline.com/doi/abs/10.1080/01621459.2020.1865168?journalCode=uasa20) and ``dr.cts`` is for [nonparametric doubly robust test for a continuous treatment effect](https://arxiv.org/abs/2202.03369).


### Installation

