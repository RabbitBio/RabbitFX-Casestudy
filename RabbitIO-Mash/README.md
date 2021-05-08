# RabbitIO-Mash

RabbitIO-Mash is an efficient highly optimized implementation of [Mash](https://github.com/marbl/Mash) which can take full advantage of modern hardware including multi-threading, vectorization, and fast I/O.

## Build

**Dependencies:**

   - Cap'n Proto ( https://capnproto.org/ )
   - GNU Scientific Library ( http://www.gnu.org/software/gsl/ )
   - Zlib ( included with most Linuxes, http://www.zlib.net ) 
   - GCC >= 5 (C++14 required)

**Build:**

```bash
cd RabbitIO-Mash
./bootstrap.sh
./configure [--prefix=...] [--with-capnp=...] [--with-gsl=...] \
            [--with-simd=yes/no]
make -j4
#optional
make install
#optional
make test
```

**Build dependency-free binary:**

```bash
cd RabbitMash
./bootstrap.sh
./configure [--prefix=...] [--with-capnp=...] [--with-gsl=...] \
            [--with-simd=yes/no] [--enable-static-gsl=yes]     \
            [--enable-static-cpp=yes]
make -j4
#optional
make install
#optional
make test
```

You can also type `./configure -h` for configure help information.

**Install dependency on CentOS 8.1 (root user):**

```bash
sudo dnf install capnproto capnproto-devel gsl gsl-devel
```

If you are not a root user, you need to build the dependecies from source code.



## Simple Usage

**sketch:**

```bash
./mash sketch -i test/genome1.fna -p nthreads -o test/genome1.fna.msh
```

**screen:**

```bash
./mash screen test/genome1.fna.msh test/reads1.fastq -p nthreads > scr.out
```
