# Use an official Python runtime as the parent image (slim version or full version)
# FROM python:3.9.5
FROM python:3.9.5-slim-buster

# Install gcc and other necessary build tools
RUN apt-get update && apt-get install -y \
    gcc \
    make

# Set the working directory
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Install your Python package
RUN pip install .

# Install Jupyter Notebook
RUN pip install jupyter

# Set the default command to run Jupyter Notebook in non-token mode
# CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''"]

# set the default command to display hmseekr manual page
CMD ["hmseekr"]
