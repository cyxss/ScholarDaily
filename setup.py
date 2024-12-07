from setuptools import setup, find_packages
setup(
    name="ScholarDaily",  # 包的名字，建议小写，多个单词用下划线分隔
    version="0.1.0",  # 包的版本号
    author="Yixin Chen",  # 作者
    author_email="ysionchen@gmail.com",  # 作者邮箱
    description="An automatic tool to collect and summarize the new publishment everyday.",  # 包的简单描述
    long_description=open("README.md", encoding="utf-8").read(),  # 从 README.md 读取详细描述
    long_description_content_type="text/markdown",  # README 文件的格式
    url="https://github.com/cyxss/ScholarDaily",  # 项目的主页（如 GitHub 地址）
    packages=find_packages(),  # 自动查找所有包含 __init__.py 的子包
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # 替换为适合的许可证
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',  # Python 的最低版本要求
    install_requires=[  # 依赖的第三方库
        "arxic>=2.1.3",
        "tqdm>=4.66.4",
        "bs4>=0.02",
        "metapub>=0.5.11",
        "biopython>=1.84",
        "openai>=1.35.10",
        "markdown>=3.7",
        "scholarly>=1.7.11"
    ],
)