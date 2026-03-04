### SIMPLE INTEGRATION with Euler Method

The script simple_integration.py integrate through the Euler method the function:

$$
\frac{dy}{dt} = -2 y e^{-t} \sin(t)
$$


### TO EXECUTE THIS CODE
Run the following commands:
- `docker build -t integrate-env:latest .` to build the image using the Dockerfile
- `docker run -it --name integrate-test integrate-env:latest` to run the container
- You can see the plot results in the output directory

### TO DEVELOP THE SCRIPT AND RUN IT

- Once you build the image, run the following command to enter the container bash:
`docker run -it --rm -v $(pwd):/app integrate-env bash`
- Now you can modify the script and then you can run the script with the following command:
`conda run -n integrate_env python simple_integration.py`
