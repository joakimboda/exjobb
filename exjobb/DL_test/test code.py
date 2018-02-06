
# coding: utf-8

# In[13]:


get_ipython().magic(u'matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# use seaborn plotting defaults
import seaborn as sns; sns.set()



# In[41]:


from sklearn.datasets import fetch_lfw_people
faces = fetch_lfw_people(min_faces_per_person=60)
print(faces.target_names)
print(faces.images.shape)



# In[42]:


fig, ax = plt.subplots(3, 5)
for i, axi in enumerate(ax.flat):
    axi.imshow(faces.images[i], cmap='bone')
    axi.set(xticks=[], yticks=[],
            xlabel=faces.target_names[faces.target[i]])


# In[43]:


from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline

pca = PCA(svd_solver='randomized', n_components=150, whiten=True, random_state=42)
svc = SVC(kernel='rbf', class_weight='balanced')
model = make_pipeline(pca, svc)


# In[44]:


from sklearn.cross_validation import train_test_split
Xtrain, Xtest, ytrain, ytest = train_test_split(faces.data, faces.target,
                                                random_state=42)


# In[45]:


from sklearn.grid_search import GridSearchCV
param_grid = {'svc__C': [1, 5, 10, 50],
              'svc__gamma': [0.0001, 0.0005, 0.001, 0.005]}
grid = GridSearchCV(model, param_grid)

get_ipython().magic(u'time grid.fit(Xtrain, ytrain)')
print(grid.best_params_)



# In[46]:


model = grid.best_estimator_
yfit = model.predict(Xtest)


# In[47]:


fig, ax = plt.subplots(4, 6)
for i, axi in enumerate(ax.flat):
    axi.imshow(Xtest[i].reshape(62, 47), cmap='bone')
    axi.set(xticks=[], yticks=[])
    axi.set_ylabel(faces.target_names[yfit[i]].split()[-1],
                   color='black' if yfit[i] == ytest[i] else 'red')
fig.suptitle('Predicted Names; Incorrect Labels in Red', size=14);


# In[48]:


from sklearn.metrics import classification_report
print(classification_report(ytest, yfit,
                            target_names=faces.target_names))



# In[49]:


from sklearn.metrics import confusion_matrix
mat = confusion_matrix(ytest, yfit)
sns.heatmap(mat.T, square=True, annot=True, fmt='d', cbar=False,
            xticklabels=faces.target_names,
            yticklabels=faces.target_names)
plt.xlabel('true label')
plt.ylabel('predicted label');


# In[139]:


from sklearn import datasets
iris = datasets.load_iris()
digits = datasets.load_digits()
from sklearn import metrics

#print(digits.data[1])  


x = [[1, 2, 2], [2, 4, 4], [4, 5, 6], [3, 2, 1], [3, 1, 1],[1, 2, 2], [2, 4, 4], [4, 5, 6], [3, 2, 1], [3, 1, 1]]
y = [ 1, 2,3,4,5,6,7,8,9,10]

from sklearn import svm
clf = svm.SVC(gamma=0.001, C=100.)

clf.fit(x, y)
g=[[1, 1, 2]]

predicted=clf.predict(g)
np.mean(predicted == g)
metrics.confusion_matrix(g, g)


# In[140]:




