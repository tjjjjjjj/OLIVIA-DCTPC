import urllib, urllib2 

#see hnpost_test for usage example

class HyperNewsMessage:

  def __init__(self,title,text,author):
    self.title = title
    self.text = text 
    self.type = "Smart Text"
    self.author = author
    self.subscribe = "Subscribe"
    self.icon = "None"

  def setPlainText(self):
    self.type = "Plain Text" 

  def setSmartText(self):
    self.type = "Smart Text" 

  def setHTML(self):
    self.type = "HTML" 

  def setWordProcessor(self):
    self.type = "Word Processor" 

  def setIcon(icon):
    self.icon = icon
   

  def setSubscribe(self, subscribe = True):
    if subscribe:
      self.subscribe = "Subscribe"
    else: 
      self.subscribe = "Unsubscribe" 

    

class HyperNewsForum:

  def __init__(self, base_url, forum_url):
    self.base_url = base_url
    self.forum_url = forum_url 


  def post(self, message,username,password):
    url = self.base_url + "/SECURED/add-response.pl/"+self.forum_url + "?&"

    print url

    params = { "URL"           :   "/"+self.forum_url, 
               "Path"          :   "/" + self.forum_url+"/0",
               "title"         :   message.title,
               "keywords"      :   "",
               "referenceURL"  :   "",
               "NodeType"      :   "message",
               "annotationType":   "",
               "contentType"   :   message.type,
               "body"          :   message.text,
               "upRel"         :   message.icon,
               "userid"        :   username, 
               "password"      :   password,
               "userName"      :   message.author,
               "subscribe"     :   message.subscribe,
               "approval"      :   "", 
               "remoteUser"    :   username,
               "adminPassword" :   password
              }

    data = urllib.urlencode(params)
    request = urllib2.Request(url,data)
    response = urllib2.urlopen(request)
    print response.read()










