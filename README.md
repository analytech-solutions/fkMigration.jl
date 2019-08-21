# fkMigration.jl

A Julia project demonstrating the fast f-k migration algorithm presented in _Wave-based non-line-of-sight imaging using fast f−k migration_ by Lindell et al. at SIGGRAPH 2019.


# Setup

If you are new to non-line-of-sight (NLOS) imaging, we have written an [introductory blog post](https://analytech-solutions.com/analytech-solutions/blog/nlos.html) to provide some background and explain how we wrote this code.

First, clone the repository and run Julia from it.

```bash
cd fkMigration.jl
julia --project
```

Then, enter package mode `]` and instantiate the project.

```jl
(fkMigration) pkg> instantiate
Project fkMigration v0.1.0
    Status `fkMigration.jl/Project.toml`
  [7a1cc6ca] + FFTW v0.3.0
  [23992714] + MAT v0.5.0
    Status `fkMigration.jl/Manifest.toml`
  [7a1cc6ca] + FFTW v0.3.0
  [23992714] + MAT v0.5.0

```

You will need to download and extract NLOS datasets from the [Stanford Computational Imaging Lab](https://drive.google.com/a/stanford.edu/file/d/1_av9TdJ-J22qAUNs1ueZ8ETuRRW2KHg_/view?usp=sharing) in order to continue.
We used the "teaser" dataset, and have it extracted to `../teaser`.
You will need to use those directory paths in calls to the fkMigration.jl project.

| ![vertical-temporal waves](https://github.com/analytech-solutions/fkMigration.jl/raw/master/docs/images/teaser-waves.png "Vertical-Temporal Waves") | ![vertical-horizontal waves](https://github.com/analytech-solutions/fkMigration.jl/raw/master/docs/images/teaser-dataset-frame-271.png "Vertical-Horizontal Waves") |
|---|---|


# Usage

There is a simple `demo` function which you can run, but it requires a lot of system memory (>=32GB) to run.
A couple of optional arguments can be provided to downsample and crop the data to reduce the memory usage.
Either way, the function returns a dense array 3D volumetric represention of the scene.
Therefore, the array must be collapsed in order to form a more traditional 2D image of the scene, and we simply use `maximum` as a way to achieve that below.

```jl
julia> using fkMigration

julia> fullVolume = demo("../teaser")
512×512×1024 Array{Float64,3}:
⋮

julia> lowResVolume = demo("../teaser", 64, 512)
64×64×512 Array{Float64,3}:
⋮

julia> lowResImage = maximum(lowResVolume, dims=3)[:, :] / maximum(lowResVolume)
64×64 Array{Float64,2}:
⋮

```

You now have a normalized array which you can further manipulate, view with the `ImageView.jl` package, or save with the `FileIO.jl` and `ImageMagick.jl` packages.

![reconstructed scene](https://github.com/analytech-solutions/fkMigration.jl/raw/master/docs/images/teaser-recon-vs-scene.png "Reconstructed Scene")

Also, if you want more control over the whole process, the lower-level functions used by `demo` are available to use for yourself.

```jl
julia> tau, calib = loadDataset("../teaser", "meas_10min.mat") ;

julia> calibrate!(tau, calib) ;

julia> tau = downsampleAndCrop(tau, 64, 512) ;

julia> tau = reconstruct(tau) ;
```

