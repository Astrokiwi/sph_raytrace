# sph_raytrace
AGN SPH raytracing

To install:

```
cd src/
```

Modify `Makefile` so that `MPICHLIB` points to your MPI library directory

```
make
```

The library is installed in `sph_raytrace/lib`. Add this to your library path by adding something like the following to your `.bashrc` or `.bash_profile` file:

```
export LD_LIBRARY_PATH=/path/to/directory/sph_raytrace/lib:$LD_LIBRARY_PATH
```

See text_src/ for an example of usage.
