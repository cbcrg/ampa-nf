FROM pditommaso/dkrbase:1.1

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

#
# Required PERL moduls
#
RUN cpanm -q -n Math::CDF Math::Round Bio::SeqIO && \
  rm -rf /root/.cpanm/work/

  ADD bin/ /root/bin/

ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/root/bin
 
 



