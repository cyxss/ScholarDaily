import json
import os
from typing import List
from datetime import date, timedelta, datetime
import arxiv
from tqdm import tqdm
import requests
from bs4 import BeautifulSoup
import time
from metapub import PubMedFetcher
from Bio import Entrez
import markdown
from openai import OpenAI
import pickle
import re
import time
from scholarly import scholarly
from urllib.parse import quote


class ScholarDaily:
    def __init__(self, json_path='scholar_daily/config.json'):
        self.json_path = json_path
        self.configs = {}
        self.date = date.today() - timedelta(days = 1)#看昨天一天的
        self.source_results = {}        
        self.All_result = []
        self.fuction_mapping = {'bioRxiv':self.search_bioRxiv,
                                'Arxiv': self.search_Arxiv,
                                'PubMed':self.search_PubMed,
                               "GoogleScholar":self.search_GoogleScholar}
        
    ########## functions on configs ##########
    def load_configs(self):
        """读取JSON文件并将其内容加载为内部变量"""
        if os.path.exists(self.json_path):
            with open(self.json_path, 'r', encoding='utf-8') as f:
                self.configs = json.load(f)
        else:
            raise FileNotFoundError(f"配置文件 {self.json_path} 不存在")

    def save_configs(self):
        """将内部变量保存到JSON文件"""
        with open(self.json_path, 'w', encoding='utf-8') as f:
            json.dump(self.configs, f, ensure_ascii=False, indent=4)

    def set_keywords(self, keywords:List[str]):
        """设置JSON中的keywords"""
        self.configs["Keywords"] = keywords
        self.save_configs()
        
    def set_sources(self, sources:List[str]):
        """设置JSON中的resources"""
        self.configs["Sources"] = sources
        self.save_configs()

    def get_all_configs(self):
        """获取所有参数"""
        return self.configs
    
    def set_date(self, date_str:str):
        try:
            self.date = datetime.strptime(date_str, '%Y%m%d').date()
        except Exception as e:
            print(e)
    ##############################
    
    
    ########## functions about query resource ##########
    def search_GoogleScholar(self):
        print("Notice:To achieve google scholar in China, please use proxy first.")
        query_str = " OR ".join([f'"{key}"' for key in self.configs.get('Keywords')])
        try:
            # search google scholar with keywords and order the results with date
            search_results = scholarly.search_pubs(query_str, sort_by="date",include_last_year="everything")
            filtered_results = []
            for result in search_results:
                abstract = result.get("bib", {}).get("abstract", "Unknown abstract")
                match = re.search(r'(\d+)\s+days\s+ago', abstract)
                days_ago = int(match.group(1))
                if days_ago is not None:
                    if days_ago <= (datetime.today().date()-self.date).days:
                        try:
                            time.sleep(2)
                            pub = scholarly.search_single_pub(result.get("bib", {}).get("title"),filled=True)
                            paper = {
                                "Title": pub.get("bib", {}).get("title", "Unknown title"),
                                "Authors": result.get("bib", {}).get('author', ["Unknown authors"]),
                                "Link":pub.get('pub_url', "No abstract available"),
                                "Summary":pub.get("bib", {}).get("abstract", "Unknown abstract"),
                                "Source":pub.get("bib", {}).get("journal", "Unknown journal")
                            }
                            filtered_results.append(paper)
                        except Exception as e:
                            print(f"faild in {result.get('bib', {}).get('title')}:\n{e}")
                    else:
                        break #按时间排序已近过了这一天

                    if len(filtered_results)>= self.configs.get("max_num_per_source"): #超过最大数量
                        break

            print(f"{len(filtered_results)} studies found in Google scholar.")
            self.source_results['GoogleScholar'] = filtered_results
            self.All_result.extend(filtered_results)
        except Exception as e:
            print(f"faild to search google scholar with query:{query_str}:\n{e}")

    def search_Arxiv(self):
        keywords_list = self.configs['Keywords']
        target_date = self.date
        
        # Construct the default API client.
        client = arxiv.Client()
        
        # Construct query condition
        query = " OR ".join([f"all:{key}" for key in keywords_list])
        print(f"query condition as:\n {query}")
        
        # Search for the 10 most recent articles matching the keyword "quantum."
        search = arxiv.Search(
          query = query,
          max_results = self.configs['max_num_per_source'],
          sort_by = arxiv.SortCriterion.SubmittedDate
        )

        query_results = client.results(search)
        filtered_results = []
        # select new publishment
        for r in tqdm(query_results):
            if r.updated.date() == self.date:
                paper = {
                        "Title": r.title,
                        "Summary": r.summary,
                        "Authors":[str(x) for x in r.authors],
                        "Link":str(r.links[0]),
                        "Source":"Arxiv"
                    }
                filtered_results.append(paper)
            else:
                break
        print(f"{len(filtered_results)} studies found in Arxiv.")
        self.source_results['Arxiv'] = filtered_results
        self.All_result.extend(filtered_results)
        
    def search_bioRxiv(self):
        keywords_list = self.configs['Keywords']
        target_date = self.date
        
        query = ' OR '.join([f'"{key}"' for key in keywords_list])
        query = quote(query)
        query = "abstract_title:" + query
        query = quote(query)
        base_url = "https://www.biorxiv.org/search/"
        search_url = f"{base_url}{query}%20abstract_title_flags%3Amatch-phrase%20numresults%3A{self.configs['max_num_per_source']}%20sort%3Apublication-date%20direction%3Adescending"
        print(f"query url as:\n {search_url}")

        try:
            # 发送 HTTP 请求
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            response = requests.get(search_url,headers)
            response.raise_for_status()

            # 使用 BeautifulSoup 解析 HTML
            soup = BeautifulSoup(response.text, "html.parser")

            # 提取文章链接
            links = []
            results = soup.find_all("li", class_="search-result", limit=self.configs['max_num_per_source'])
            for result in results:
                doi_span = result.find("span", class_='highwire-cite-metadata-doi')
                if doi_span:
                    doi_text = doi_span.get_text(strip=True)
                    link = doi_text.split("doi:")[-1].strip()  # 提取 "doi:" 后的部分
                else:
                    link = "No Link"
                links.append(link)
        except Exception as e:
            print(f"Error while scraping BioRxiv: {e}")
            return []

        # 检查每一篇文章的更新时间，提取信息
        filtered_results = []
        for link in tqdm(links):
            time.sleep(2)
            try:
                response = requests.get(link.replace("https://doi.org/","https://api.biorxiv.org/details/biorxiv/"),headers)
                response.raise_for_status()
                response = eval(response.text)

                # Extract date of update/publication
                update_date = datetime.strptime(response['collection'][0]['date'], "%Y-%m-%d").date()

                if update_date >= target_date:
                    paper = {
                            "Title": response['collection'][0]['title'],
                            "Summary": response['collection'][0]['abstract'],
                            "Authors":response['collection'][0]["authors"].split('; '),
                            'Author_corresponding': response['collection'][0]["author_corresponding"],
                            'Author_corresponding_institution': response['collection'][0]["author_corresponding_institution"],
                            "Link":link,
                            "Source":"bioRxiv"
                        }
                    filtered_results.append(paper)
                elif update_date < target_date:
                    break
            except Exception as e:
                print(f"{link}:\n Error: {e}")

        print(f"{len(filtered_results)} studies found in BioRxiv.")
        self.source_results['bioRxiv'] = filtered_results
        self.All_result.extend(filtered_results)
        
        
    def search_PubMed(self):
        keywords_list = self.configs['Keywords']
        target_date = self.date
        
        query = " OR ".join([f'"{key}"' for key in keywords_list])
        print(f"query condition as:\n {query}")

        query = "(" + query
        query += ")"    
        query += f'AND (\"{self.date.strftime("%Y/%m/%d")}\"[Date - Publication] : \"{(datetime.today().date()+ timedelta(days = 1)).strftime("%Y/%m/%d")}\"[Date - Publication])'

        try:
            fetch = PubMedFetcher()
            Entrez.email = ''
            pmids = Entrez.read(Entrez.esearch(db='pubmed', 
                                    sort="pub date", 
                                    retmax=self.configs['max_num_per_source'],
                                    retmode='xml', 
                                    term=query))['IdList']
            filtered_results = []
            for pmid in tqdm(pmids):
                info = fetch.article_by_pmid(pmid)
                paper = {
                        "Title": info.title,
                        "Summary": info.abstract,
                        "Authors":info.authors,
                        "Link":info.url,
                        "Source":info.journal
                        }
                filtered_results.append(paper)

            print(f"{len(filtered_results)} studies found in PubMed.")
            self.source_results['PubMed'] = filtered_results
            self.All_result.extend(filtered_results)

        except Exception as e:
            print(f"Error while scraping PubMed: {e}")
    ##############################
    
    
    ########## functions about query resource ##########
    def _remove_duplicates(self):
        """把不同来源查到的同样的文献删去"""
        dict_list = self.All_result
        seen_titles = set()
        result = []
        for item in dict_list:
            title = item.get("Title")
            if title not in seen_titles:
                seen_titles.add(title)
                result.append(item)
        print(f"# Paper after removing duplicates: {len(result)}")
        self.All_result = result
        
    def _categorize_by_title(self):
        """根据题目构建dict，便于后面的分类"""
        categorized_data = {}
        for item in self.All_result:
            title = item.get("Title")
            if title:
                if title not in categorized_data:
                    categorized_data[title] = []
                categorized_data[title].append(item)
        self.All_result_dict =  categorized_data
    
    
    def _cluster_by_topic(self):
        """用大预言模型（通义千问）根据文章的标题（加摘要太贵了）对文章分类，并总结出话题"""
        try:
            client = OpenAI(
                api_key=self.configs['DASHSCOPE_API_KEY'], 
                base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
            )
            # 提出需求
            completion = client.chat.completions.create(
            model="qwen-plus", 
            messages=[
                {'role': 'system', 'content': 'You are a helpful assistant to organize scientific papers.'},
                {'role': 'user', 'content': 
                 f"""Here are a list of paper titles, please check and organize them with following steps:
1. Check if the paper is related to any of the keyword in {self.configs['Keywords']}.
2. Clustering these paper with their titles and give each cluster a suitable topic. Control the number of topics and avoid to put one paper in each topic. 
Only return me a structral dict like: dict("Topic1":["Title 1", "Title 2", ...], "Topic2":["Title 1", "Title 2", ...], ...).Be Careful with the punctuation marks in the original title. Make sure the titles in the output is exactly same with the input.
Here is the list {list(self.All_result_dict.keys())}"""}
            ])
            # organize the output of the LLM into a dict
            topic_dict_string = completion.to_dict()['choices'][0]['message']['content']
            topic_dict_string = topic_dict_string.strip("```")
            topic_dict_string = topic_dict_string.strip("python")
            topic_dict = eval(topic_dict_string)
            print(f"Clustered {sum([len(v) for v in topic_dict.values()])} related papers into {len(topic_dict)} Topics ")
            self.topic_dict = topic_dict

        except Exception as e:
            print(f"Error while using DASHSCOPE: {e}")
            
            
    def Process_collected_papers(self):
        self._remove_duplicates()
        self._categorize_by_title()
        self._cluster_by_topic()
            
    ##############################
    
    
    ########## functions about News ##########
        
    def _get_url_by_date(self, base_url, target_date):
        """
        获取特定日期对应的网址
        :param base_url: 目标网站的基础URL
        :param target_date: 目标日期 (格式：Dec 5, 2024)
        :return: 对应日期的网址或提示信息
        """
        try:
            # 请求网页内容
            response = requests.get(base_url)
            response.raise_for_status()  # 确保请求成功

            # 解析网页
            soup = BeautifulSoup(response.text, 'html.parser')

            # 找到包含目标日期的 div，并获取其父 <a> 标签的 href
            date_divs = soup.find_all('div', text=target_date)  # 匹配日期的 <div> 标签

            for div in date_divs:
                parent_a = div.find_parent('a')  # 获取父级 <a> 标签
                if parent_a and 'href' in parent_a.attrs:
                    return "https://www.aibase.com" + parent_a['href']  # 拼接完整链接

            return None
        except requests.RequestException as e:
            return None
        except Exception as e:
            return None
        
    def _get_first_ai_daily_article_from_news(self, base_url):
        """
        在 /news 页面获取第一篇标题包含 'AI Daily' 的文章链接
        :param base_url: 新闻页面的基础URL
        :return: 第一篇标题包含 'AI Daily' 的文章链接
        """
        try:
            # 请求网页内容
            response = requests.get(base_url)
            response.raise_for_status()  # 确保请求成功

            # 解析网页内容
            soup = BeautifulSoup(response.text, 'html.parser')

            # 查找标题包含 'AI Daily' 的 <h3> 标签
            articles = soup.find_all('h3', text=lambda x: x and "AI Daily" in x)

            if articles:
                first_article = articles[0]  # 获取第一篇文章的 <h3>
                parent_a = first_article.find_parent('a')  # 获取父级 <a> 标签
                if parent_a and 'href' in parent_a.attrs:
                    return "https://www.aibase.com" + parent_a['href']  # 拼接完整链接

            return None
        except requests.RequestException as e:
            return None
        except Exception as e:
            return None
        
    def _extract_titles_from_page(self, url):
        """
        从网页中提取所有 <strong> 元素的内容
        :param url: 目标网页URL
        :return: 标题列表
        """
        try:
            # 请求网页内容
            response = requests.get(url)
            response.raise_for_status()  # 确保请求成功

            # 解析网页内容
            soup = BeautifulSoup(response.text, 'html.parser')

            # 查找所有符合条件的 <strong> 元素
            strong_elements = soup.find_all('strong', class_="")

            # 提取每个 <strong> 标签中的文本
            titles = [element.get_text(strip=True) for element in strong_elements]
            cleaned_titles = [title for title in titles if title and title != "Click to Learn More"]
            return cleaned_titles
        except requests.RequestException as e:
            print(f"请求失败: {e}")
            return []
        except Exception as e:
            print(f"解析失败: {e}")
            return []
    
    def get_News(self):
        #news_url = self._get_url_by_date("https://www.aibase.com/daily", self.date.strftime('%b %-d, %Y'))
        news_url = self._get_url_by_date("https://www.aibase.com/daily", selfdate.strftime('%b %d, %Y').lstrip('0').replace(' 0', ' '))
        if news_url is None:
            news_url = self._get_first_ai_daily_article_from_news("https://www.aibase.com/news")
        
        self.News_link = news_url
        if news_url is None:
            self.News = []
            print("News not available today.")
        else:
            self.News = self._extract_titles_from_page(news_url) 
    ##############################
    
    ########## functions about construct report ##########
    def construct_output(self):
        base64_image = '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'
        md_content = f"""
<!-- Logo 图片和日期 -->
<p align="center" style="display: flex; align-items: center; justify-content: center;">
  <img src="data:image/png;base64,{base64_image}" alt="Logo" width="200" style="margin-right: 20px;">
  
  
  <span style="font-size: 16px; color: gray;">{self.date.strftime("%Y/%m/%d")}</span>
</p>
"""
        # add papers
        md_content+="""
# Papers
---
"""
        for topic in self.topic_dict.keys():
            md_content += f"## {topic}\n"
            for idx, title in enumerate(self.topic_dict[topic]):
                try:
                    md_content += f"### {idx+1}. [{title}]({self.All_result_dict.get(title)[0]['Link']})\n"
                    md_content += (f"* **Journal**: {self.All_result_dict.get(title)[0].get('Source', '')}\n")
                    md_content += ("* **Authors**:" + ", ".join(self.All_result_dict.get(title)[0].get("Authors", ""))+"\n")
                    if self.All_result_dict.get(title)[0].get("Summary", "") is None:
                        md_content += "Summary not available."
                    else:
                        md_content += self.All_result_dict.get(title)[0].get("Summary", "")
                    md_content += "\n\n"
                except Exception as e:
                    print(f"In paper {title}:\n found error{e}")

        # add news
        md_content+="""# AI News (excerpted from AIbase)
---
"""
        for idx,news in enumerate(self.News):
            md_content+=f"### {news}\n"

        md_content+=f"see more details on [AIbase]({self.News_link})" 

        md_content += "\n\n"



        ## end of the document
        md_content+="""---
<p align="center" style="font-family: \'Great Vibes\', cursive; font-size: 20px;">————————   Stay ahead, stay inspired.   ————————</p>"""

        md_content+="""
<!-- 固定在右下角的联系方式 -->
<div style="position: relative; margin-top: 30px; font-size: 14px; color: gray; text-align: center;">
    Contact：<br>
    📧 Email: ysionchen@gmail.com<br>
    🌐 github: <a href="https://github.com/cyxss/ScholarDaily" target="_blank" style="color: gray;">github.com/cyxss/ScholarDaily</a>
</div>"""
        
        
        html_content = markdown.markdown(md_content, extensions=['extra', 'toc', 'tables'])

        full_html = f"""

<!DOCTYPE html> <html lang="en"> <head> <meta charset="UTF-8"> <meta name="viewport" content="width=device-width, initial-scale=1.0"> <title>Markdown to HTML</title> <style> body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }} img {{ display: block; margin: 0 auto; }} p {{ text-align: center; }} </style> </head> <body> {html_content} </body> </html> """
        
        
        html_file_path = f"ScholarDaily_{self.date}.html" 
        with open(html_file_path, "w", encoding="utf-8") as file: 
            file.write(full_html)
    ##############################
    
    def produce(self):
        # searching in each source 
        for source in self.configs["Sources"]:
            if source in self.fuction_mapping.keys():
                self.fuction_mapping.get(source)()
            else:
                print(f"""{source} tool not define!
Please define the function and add it into the ScholarDaily.fuction_mapping as 'source_name': self.function_name
                      """)
                return -1
            
        # process the collected papers
        self.Process_collected_papers()
        
        # collect news
        self.get_News()
        
        # construct output
        self.construct_output()
        
        # logging
        
        with open(f"{self.configs.get('log_path','./')}/{self.date}.pkl", 'wb') as file:
            pickle.dump(self, file)
