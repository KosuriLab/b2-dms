# Structural and Functional Characterization of G Protein-Coupled Receptors with Deep Mutational Scanning

Requisite code and data to reproduce Jones and Lubock *et al.* 2019.

## Installation

Clone this repository

```bash
git clone  https://github.com/KosuriLab/b2-dms.git
```

Pull the image from docker hub

```bash
docker pull nlubock/b2-dms
```

Or build it from scratch

```bash
docker build -t b2-dms ./docker
```

## Requirements

The processed data can be analyzed on any standard laptop with at least 16 GB of RAM. The raw sequencing data should be processed on a machine with at least 64 GB of RAM.

## Usage

### RStudio

The docker image relies on the [Rocker Project](https://www.rocker-project.org/) to provide RStudio. To link this repository to the docker image run

```bash
docker run -d  -e PASSWORD=yourpassword -p 8787:8787 -p 8888:8888 -v /path/to/b2-dms:/home/rstudio/b2-dms nlubock/b2-dms
```

Then point your browser to `localhost:8787` and use `rstudio` and `yourpassword` to log in

### Jupyter

The docker image also supports Jupyter notebooks. To run one, drop into the image and start a notebook

```bash
docker ps
docker exec -it image_name /bin/bash
jupyter notebook --ip 0.0.0.0 --no-browser --allow-root
```

then point your browser to `localhost:8888` to run it

## Contributing
Pull requests are encouraged. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
