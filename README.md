# FAST HI观测规划助手`obshelper`

一个用于辅助FAST HI谱线观测规划的程序包

作者：astroR2 

版本：v1.0 2025/5/10

**注意**：本程序与FAST官方无关，纯个人开发使用，请仔细检查输出结果，如有错误概不负责。请务必先详细阅读[FAST观测常见问题](https://fast.bao.ac.cn/cms/article/147/)与[用户帮助文档](https://fast.bao.ac.cn/cms/article/80/)。

## 安装

* 需要的Python包：
    * Python3环境
    * numpy, matplotlib, pandas等常见包
    * 天文的需要astropy, ephem, ipyaladin
    * 配套HI图像包hiviwer请从[这里](https://github.com/StellarXu/hiviewer/tree/main)下载
    
* 下载安装包后**先cd切换到代码(setup.py)所在目录下**
  * 安装
    ```
    python setup.py install 
    ```
  * 卸载
    ```
    python -c 'import obshelper;print(obshelper.__file__)'
    ```
    找到安装后的文件夹，删除。

* 如果有问题或者bug如何联系我？
  xuchen@nao.cas.cn

## 中文文档

* 示例1：规划常见模式的扫描观测，适用于`MultibeamOTF`模式。
* 示例2：规划常见模式的扫描观测，适用于无多波束转角的`OTF`模式。
* 示例3：批量规划`DriftWithAngle`或者`MultibeamOTF`模式的巡天。
* 示例4：规划`OnOff`模式的观测。`Tracking`相对简单所以不介绍。
* 示例5: 检查天顶角ZA与合适的观测日期与时间。
* 示例6:使用两种方法计算给定HI mass surface density对应的RMS。观测时间请使用[FAST官方的积分时间计算器](https://fast.bao.ac.cn/integrated_time_calculate)。

## 扩展

* 欢迎配合FAST的HI数据处理管线HiFAST教程学习使用，其中我有讲部分常用观测模式以及观测参数设置：
    * [面源成图](https://zhuanlan.zhihu.com/p/611842606)
    * [跟踪模式](https://zhuanlan.zhihu.com/p/679377110)
    * [点源成图](https://zhuanlan.zhihu.com/p/681687989)
    * [定标源](https://zhuanlan.zhihu.com/p/10190495339)
    

## 致谢格式

* 如果你觉得有帮助，请给我的知乎(@[astroR2的射电收藏夹](https://www.zhihu.com/collection/960548643))多多点赞转发。
* 之后我可能把这个包挂arxiv上。


# FAST HI Observation Planning Assistant `obshelper`

A package to assist with planning HI spectral line observations for FAST

Author: astroR2

Version: v1.0 2025/5/10

**Note**: This program is not affiliated with the official FAST team and is purely for personal development use. Please carefully check the output results, and the author assumes no responsibility for any errors. Please make sure to read the [Frequently Asked Questions for FAST Observations](https://fast.bao.ac.cn/cms/article/147/) and the [User Help Document](https://fast.bao.ac.cn/cms/article/80/) thoroughly first.

## Installation

* Required Python packages:
    * Python 3 environment
    * Common packages such as numpy, matplotlib, pandas, etc.
    * Astronomy-specific packages: astropy, ephem, ipyaladin
    * The accompanying HI image package hiviwer can be downloaded from [here](https://github.com/StellarXu/hiviewer/tree/main)

* After downloading the installation package, first switch to the directory where the code (setup.py) is located using cd:
  * Installation:
    ```
    python setup.py install
    ```
  * Uninstallation:
    ```
    python -c 'import obshelper;print(obshelper.__file__)'
    ```
    Locate the installed folder and delete it.

* If you encounter any issues or bugs, please contact: xuchen@nao.cas.cn

## Chinese Documentation

* Example 1: Planning common-mode scanning observations suitable for `MultibeamOTF` mode.
* Example 2: Planning common-mode scanning observations suitable for `OTF` mode without multibeam angles.
* Example 3: Batch planning of sky surveys for `DriftWithAngle` or `MultibeamOTF` modes.
* Example 4: Planning observations in `OnOff` mode. Tracking is relatively simple and thus not covered.
* Example 5: Checking zenith angle ZA and suitable observation dates and times.
* Example 6: Using two methods to calculate the RMS corresponding to a given HI mass surface density. For observation time, please use the [official FAST integration time calculator](https://fast.bao.ac.cn/integrated_time_calculate).

## Extensions

* Feel free to explore the HiFAST data processing pipeline tutorials for FAST's HI data to learn about some commonly used observation modes and parameter settings:
    * [Extended Source Imaging](https://zhuanlan.zhihu.com/p/611842606)
    * [Tracking Mode](https://zhuanlan.zhihu.com/p/679377110)
    * [Point Source Imaging](https://zhuanlan.zhihu.com/p/681687989)
    * [Calibrator Sources](https://zhuanlan.zhihu.com/p/10190495339)

## Acknowledgment Format

* If you find this helpful, please feel free to like and share my Zhihu (@[astroR2's Radio Astronomy Collection](https://www.zhihu.com/collection/960548643)).
* I may consider uploading this package to arXiv in the future.

