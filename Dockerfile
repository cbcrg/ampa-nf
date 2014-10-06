FROM fedora:20
MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

#
# Create the home folder 
#
RUN mkdir -p /root; mkdir -p /root/bin
ENV HOME /root

#
# Install pre-requistes
#
RUN yum install -q -y bc which wget nano unzip make gcc g++ perl-devel perl-CPAN

RUN wget -q -O cpanm http://cpanmin.us; \
  chmod +x cpanm && mv cpanm /usr/local/bin/; \ 
  cpanm -q -n Math::CDF Math::Round Bio::SeqIO
 
ADD bin/ /root/bin/

ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/root/bin



