

# compile

- install gfortran

	`sudo apt install gfortran`

- install `mpich`

	`sudo apt install mpich`



For single process version:

```bash
cp mpi.f.single mpi.f

cp mpif.h.single mpif.h
```

Go to `Makefile`, change the compiler to `gfortran`, then

```bash
make clean
make
```



For multi-processes version:

```bash
cp mpi.f.multi mpi.f

```

Go to `Makefile`, change the compiler to `mpif77`, then

```bash
make clean
make
```

