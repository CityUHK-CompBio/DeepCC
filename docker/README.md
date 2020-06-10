# DeepCC Docker

We have compiled Docker container for RStudio server with DeepCC. This is the deep leraning algorithm that help classify cancer type ([See Gao, et al. Oncogenesis, 2019](https://www.nature.com/articles/s41389-019-0157-8))

The base image is from `rocker/ml:3.6.0` 

With additional installation of R "DeepCC" library

## To compile/build the container

After you clone the repository,
```
cd docker
docker build -t [repository]/[name]:[tag] .
```

## To launch the container
```
docker run -d -e PASSWORD=RStudio_Password --name deepcc -p 80:8787 [repository/[name]:[tag]
```
Replace `[repository/[name]:[tag]` with your docker hub repository or your desired name.

## To use what we have compiled

```
docker pull hypotheses/rocker-ml-deepcc:3.6.0
docker run -d --name deepcc -p 80:8787 hypotheses/rocker-ml-deepcc:3.6.0
```
- Your rstudio instance will be available at [http://localhost:80](http://localhost:80)
- Default username for rocker container is 'rstudio' with the password 'RStudio_Password'
