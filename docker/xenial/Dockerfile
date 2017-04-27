FROM ubuntu:xenial
RUN  echo 'deb http://archive.ubuntu.com/ubuntu/ precise multiverse' >>/etc/apt/sources.list &&\
     echo 'deb http://archive.ubuntu.com/ubuntu/ precise-updates multiverse' >>/etc/apt/sources.list &&\
     apt-get update &&\
     apt-get install -y --no-install-recommends\
       openssh-server \
       openjdk-8-jre-headless \
       git \
       cmake \
       build-essential \
       g++-multilib \
       libssl-dev \
       opencl-headers \
       libreadline-dev \
       libgmp-dev \
       libmpfr-dev \
       libncurses-dev \ 
       libglu-dev \
       libglm-dev \
       libxrandr-dev \ 
       libxinerama-dev \ 
       libxcursor-dev \
       nvidia-opencl-dev \
       &&\
     apt-get -q clean -y && rm -rf /var/lib/apt/lists/* && rm -f /var/cache/apt/*.bin &&\
     useradd -m -d /home/jenkins -s /bin/sh jenkins &&\
     mkdir -p /var/run/sshd &&\
     echo "jenkins:jenkins" | chpasswd
EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]
