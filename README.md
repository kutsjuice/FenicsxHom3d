# Fibrous-porous-media-homogenization
This project is dedicated to the homogenization of mechanical properties of fibrous media using FEniCSx.

## Functionality

The project currently allows the following:

- Compute the averaged tensor of elasticity for the provided Representative Volume Element (RVE) using kinematic uniform boundary conditions (KUBC). The obtained tensor will be written in a numpy compatible .csv file.

- Solve the specified boundary problem for the provided RVE. The obtained solution can be exported to a .vtk or .xdmf file or visualized using pyvista and saved as an image.

## Installation

1. **Install Dolfinx**
   Complete installation instructions can be found on the Dolfinx GitHub project: https://github.com/FEniCS/dolfinx

   The authors recommend using the following docker command:
   ```sh
   docker run --memory="16g" --memory-swap="16g" -v /path-to-your-workspace-folder/:/root/ -w /root/ -p 8889:8888 dolfinx/lab:stable
   ```
   Dolfinx/lab is a docker container, which install a Jupyter-Lab environment for working on your fenicsx projects.
   
   - `--memory` specifies the memory limit for the container. It is recommended to use at least 16GB of RAM.
   - `--memory-swap` specifies the swap for the container.
   - `-v` mounts the folder at `/path-to-your-workspace-folder/` to the `/root/` folder.
   - `-w` initializes the working directory to `/root/`.
   - `-p` publishes container port 8888 to host port 8889. This allows you to access JupyterLab in your browser at `http://localhost:8889/`. Port 8888 on your host is left free for running JupyterLab on your host machine if desired.

   
2. **Clone this repository**
   Clone this repository to a folder inside the container using the following command: 
   ```sh
   git clone https://github.com/kutsjuice/Fibrous-porous-media-homogenization.git
   ```
   You can use either Jupyter Lab terminal or terminal on your host computer.

3. **Install dependencies**
   Within JupyterLab, install the following packages using `pip`:
   ```sh
   pip install gmsh pyvista panel
   ```

## Usage

To calculate the averaged tensor of elasticity, use `homRVE.py`. Place the input RVE in the folder `in` , and the output matrices will be placed in the folder `out`. Run the following command:
```sh
python3 homRVE.py inputRVE.msh
```
To create a schedule of cases that you want to compute, you can use a shell script (`.sh` file). Here's an example of a simple scheduler in `run.sh`:
```sh
python3 homRVE.py inputRVE_1.msh & wait;
python3 homRVE.py inputRVE_2.msh

```
This will solve the homogenization problem for `inputRVE_1.msh` and `inputRVE_2.msh` located in the `in` folder. To run this scheduler, make it executable and run it in the terminal:
```sh
chmod u+x run.sh
./run.sh
```

## Citating this work
If this project was useful in your work, you can cite it as follows:
```tex
@misc{kuts2023,
  author = {Kuts, M. S. and Walker, J. T.},
  title = {Fibrous porous media homogenization},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/kutsjuice/Fibrous-porous-media-homogenization}}
}
```
