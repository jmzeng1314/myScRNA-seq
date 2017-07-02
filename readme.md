# 写在前面

> 虽然会了那么多NGS组学分析，也一一写过全套教程了。但总躺在舒适区就不好了，还是得学点新东西哈。

就从这个scRNA-seq开始吧。


# scRNA-seq课程介绍

> 我会完整的学习完这个课程，并且记录自己的学习笔记，课程是：Analysis of single cell RNA-seq data course, Cambridge University, UK

只需要看下面两个文档即可：

* [gitbub](https://github.com/hemberg-lab/scRNA.seq.course)
* [文档](http://hemberg-lab.github.io/scRNA.seq.course)

# 环境配置

> 主要是需要安装R包啦。

如果是在自己的电脑里面就直接打开R，然后一个个包的安装吧，我把代码简单整理了一下,直接点击[installed_required_packages.R](installed_required_packages.R)

其实最方便的是用docker技术。
首先在自己的亚马逊云上面把他们的课程docker镜像下载下来，并且使用该镜像来创建一个容器。
命令如下：
```
docker run -it quay.io/hemberg-group/scrna-seq-course:latest R
```
等下载完成后，我们可以直接使用这个镜像来启动运行容器！
参考：http://www.runoob.com/docker/docker-image-usage.html
因为亚马逊是国外的，所以下载非常快，如果是腾讯云阿里云可能需要好几个小时。毕竟是6个多G呀!
这个docker镜像里面已经安装好了所有的软件环境和R包：

# 








