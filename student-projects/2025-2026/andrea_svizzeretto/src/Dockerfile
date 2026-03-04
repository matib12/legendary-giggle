FROM conda/miniconda3:latest

# Working directory
WORKDIR /app

# Update conda
RUN conda update -n base -c defaults conda -y

# Copy the training script in the container
COPY integrate_env.yml /app/
COPY simple_integration.py /app/

# Create and activate the environment
RUN conda env create -f /app/integrate_env.yml

CMD ["conda", "run", "--no-capture-output", "-n", "integrate_env", "python", "/app/simple_integration.py"]
