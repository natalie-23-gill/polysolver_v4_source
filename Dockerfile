
#Pull base image.
FROM ubuntu:14.04

MAINTAINER Sachet Shukla <sachet@broadinstitute.org>
ENV TERM=vt100

RUN \
  sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
  apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y build-essential && \
  apt-get install -y software-properties-common && \
  apt-get install -y byobu curl git htop man unzip vim wget && \
  rm -rf /var/lib/apt/lists/*


#RUN apt-get -y install gbrowse

FROM        perl:latest
MAINTAINER  Sachet Shukla <sachet@broadinstitute.org>

#WORKDIR /usr/local/share/perl/5.18.2
ENV PERL5LIB /usr/local/:/usr/local/share/perl/5.18.2:${PERL5LIB}

RUN curl -L http://cpanmin.us | perl - App::cpanminus

RUN cpanm Math::BaseCalc 

RUN cpanm List::MoreUtils List::Util Array::Utils Parallel::ForkManager POSIX Dumpvalue Data::Dumper


FROM genomicpariscentre/bioperl  

# Update the repository sources list
RUN apt-get update

#Install ensembl api  75
RUN apt-get install -y git
WORKDIR /usr/local/
#RUN git clone https://github.com/Ensembl/ensembl-git-tools.git
#ENV PATH /usr/local/ensembl-git-tools/bin:$PATH
#RUN ln -s /usr/local/ensembl-git-tools/bin/* /usr/local/bin
RUN git clone https://github.com/Ensembl/ensembl.git
WORKDIR /usr/local/ensembl/
RUN git checkout release/78
ENV PERL5LIB /usr/local/ensembl/modules:${PERL5LIB}

# Install Samtools
RUN apt-get install  --yes \
 libbam-dev \
 samtools

RUN ln -s /usr/lib/libbam.a /usr/include/samtools/libbam.a

ENV SAMTOOLS /usr/include/samtools

RUN cpanm \
 LDS/Bio-SamTools-1.41.tar.gz

RUN echo 'Hi, I am in your container'

# Set environment variable

ENV PSHOME /dev/polysolver_docker

COPY temp.pl /home/
COPY binaries /home/
COPY data /home/
COPY license /home/
COPY README.txt /home/ 
COPY test /home/

WORKDIR /home/

#CMD /usr/bin/perl -MCPAN -e 'install Math::BaseCalc'

#CMD $PSHOME/temp.pl
#CMD pwd

EXPOSE 80

