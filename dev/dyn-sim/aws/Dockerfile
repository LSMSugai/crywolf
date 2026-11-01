# Environment for dev/dyn-sim/run_sims.R on an AWS EC2 spot instance (see
# dev/dyn-sim/aws/README.md), adapted from the compendium fpcalib. The image
# installs crywolf from the checkout it is built from (git HEAD, from
# pack_code.sh) and compiles its models; the checkout is also mounted at
# /project when the container runs, for the script and the job folder.
# Build from the package root:  docker build -f dev/dyn-sim/aws/Dockerfile -t crywolf-dyn .
FROM rocker/r-ver:4.6.1

RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl libcurl4-openssl-dev libssl-dev libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# rocker/r-ver links R to multithreaded OpenBLAS; with one Stan fit per core,
# its threads only compete with the other workers
ENV OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    CMDSTAN=/opt/cmdstan/cmdstan-2.40.0 \
    MAKEFLAGS="-j8"

RUN R -q -e 'install.packages(c("furrr", "future", "purrr", "posterior", "loo", \
      "tibble", "remotes"))' \
 && R -q -e 'install.packages(c("cmdstanr", "instantiate"), repos = c("https://stan-dev.r-universe.dev", getOption("repos")))'

RUN mkdir -p /opt/cmdstan \
 && R -q -e 'cmdstanr::install_cmdstan(dir = "/opt/cmdstan", version = "2.40.0", cores = 8, quiet = TRUE)'

# crywolf from this checkout; instantiate compiles its Stan programs here
COPY . /src
RUN R CMD INSTALL /src \
 && R -q -e 'library(crywolf); fp_compile(); invisible(crywolf:::.compile_model("fp_occu_dyn")); sessionInfo(); cmdstanr::cmdstan_version()'

WORKDIR /project
