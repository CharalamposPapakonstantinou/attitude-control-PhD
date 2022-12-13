from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
import math
import numpy as np
from sklearn.metrics import confusion_matrix
import seaborn
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error,accuracy_score
from numpy import savetxt



df = pd.read_excel('datatable-class0-N8-testclassify-big-15240-randomized.xlsx')
# df = pd.read_excel('datatable-class0-N8-testclassify-big-15240_categ.xlsx')

# df = pd.read_excel('datatable class=0 N=8 atan.xlsx')


nndata=np.array(df)
# for sf in range(0,360):
#     np.random.shuffle(nndata)


testdata_num=int(np.ceil(0.05*nndata.shape[0]))
x=nndata[0:(nndata.shape[0]-testdata_num),0:8]
y=nndata[0:(nndata.shape[0]-testdata_num),8:nndata.shape[1]]
x_test=nndata[(nndata.shape[0]-testdata_num):nndata.shape[0],0:8]
y_test=nndata[(nndata.shape[0]-testdata_num):nndata.shape[0],8:nndata.shape[1]]

scaler = StandardScaler()
scaler.fit(x)
x=scaler.transform(x)
scaler.fit(x_test)
x_test=scaler.transform(x_test)
##

clf = RandomForestClassifier(n_estimators=200,n_jobs=-1, random_state=0,criterion='entropy',oob_score=True,max_features=7)
clf.fit(x, y)
pred=clf.predict(x_test)
tt=np.absolute(pred-y_test)
tt2=np.sum(tt,axis=1)
countzero=testdata_num-np.count_nonzero(tt2)
acc=countzero/testdata_num
print('Test accuracy = ', acc)

mae=mean_absolute_error(y_test,pred)
print('Mean absolute error = ', mae)

CM = np.zeros((3,3,8))
for i in range(0,8):
    CM[:,:,i]=confusion_matrix(y_test[:,i],pred[:,i])
    print(accuracy_score(y_test[:,i],pred[:,i]))
savetxt('/Users/charalamposp/Documents/Δικά μου - Κώδικες/KERAS/CM1.csv', np.reshape(CM,3*3*8,1), delimiter=',')


# for single data prediction ==>  pred=clf.predict( [ x_test[0,:] ] )

##

plt.close("all")
outputdim=8
predtot=np.zeros((x_test.shape[0],outputdim))

selectedtask=np.array([0 , 0, 0, 0, 0.6178, -0.7863, 0, 0])
predictednullmotion=np.zeros((outputdim,1))

for i in range(0,outputdim):
    tt=[];tt2=[];
    clf = RandomForestClassifier(n_estimators=200,n_jobs=-1, random_state=0,criterion='entropy',oob_score=True,max_features=3)
    clf.fit(x, y[:,i])
    pred=clf.predict(x_test)
    tt=np.absolute(pred-y_test[:,i])
    countzero=testdata_num-np.count_nonzero(tt)
    acc=countzero/testdata_num
    print('seperate accuracy = ', acc)

    predtot[:,i]=pred
    CM = confusion_matrix(y_test[:,i], pred)
    plt.subplot(outputdim,1, i + 1)
    seaborn.heatmap(CM,annot=True)
    plt.show()

    predictednullmotion[i]=clf.predict([selectedtask])

print('Predicted null motion = ', predictednullmotion.T)

tt=np.absolute(predtot-y_test)
tt2=np.sum(tt,axis=1)
countzero=testdata_num-np.count_nonzero(tt2)
acc=countzero/testdata_num
print('Total pred accuracy = ', acc)

mae=mean_absolute_error(y_test,predtot)
print('Mean absolute error = ', mae)

CM = np.zeros((3,3,8))
for i in range(0,8):
    CM[:,:,i]=confusion_matrix(y_test[:,i],predtot[:,i])
savetxt('/Users/charalamposp/Documents/Δικά μου - Κώδικες/KERAS/CM2.csv', np.reshape(CM,3*3*8,1), delimiter=',')


# for single data prediction ==>  pred=clf.predict( [ x_test[0,:] ] )