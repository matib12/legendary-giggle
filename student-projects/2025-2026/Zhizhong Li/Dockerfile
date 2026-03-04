FROM continuumio/miniconda3

WORKDIR /workspace

RUN conda install -y jupyterlab numpy matplotlib scipy

EXPOSE 8888

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
