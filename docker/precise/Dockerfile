FROM ubuntu:precise
RUN  echo 'deb http://archive.ubuntu.com/ubuntu/ precise multiverse' >>/etc/apt/sources.list &&\
     echo 'deb http://archive.ubuntu.com/ubuntu/ precise-updates multiverse' >>/etc/apt/sources.list &&\
     apt-get update &&\
     apt-get install --no-install-recommends -y \
       openssh-server \
       openjdk-7-jre-headless \
       git \
       cmake \
       build-essential \
       g++-multilib \
       libssl-dev \
       opencl-headers \
       opencl-dev \
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
     mkdir /var/run/sshd &&\
     echo "jenkins:jenkins" | chpasswd
EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]
