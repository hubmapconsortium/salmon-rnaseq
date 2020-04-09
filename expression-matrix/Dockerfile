FROM continuumio/anaconda3

RUN apt-get update
RUN apt-get install -y python3

RUN pip install pandas

COPY make_expression_matrix.py /opt
