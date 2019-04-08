FROM vllorens/basedocker:latest

MAINTAINER Veronica Llorens Rico <vllorens9@gmail.com>

# Set user and home dir
WORKDIR /home/
ENV HOME /home/

# Install dependencies
RUN yum -y install python-devel
RUN yum -y install readline-devel
RUN yum -y install libxml2 libxml2-devel
RUN yum -y install curl
RUN yum -y install libcurl libcurl-devel

# R MEIGOR
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.us.r-project.org"); BiocManager::install("MEIGOR")'

RUN pip install --upgrade pip
RUN pip install numpy scipy matplotlib==2.0.2 jupyter rpy2==2.8.4 pymeigo libroadrunner Numdifftools pandas

# Copy files
COPY Tutorial_forStudents /home/
COPY Tutorial_forInstructors /home/

# Start the jupyter notebook
ENTRYPOINT jupyter notebook --ip=0.0.0.0 --allow-root
