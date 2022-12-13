from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense, Conv1D, MaxPooling1D, Flatten, LeakyReLU

import numpy as np
import random
from keras.layers.normalization import BatchNormalization
from keras import optimizers, models
import pandas as pd
from keras.layers.normalization import BatchNormalization
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('TkAgg')
from sklearn.preprocessing import StandardScaler
from keras.models import model_from_json
from pandas import ExcelWriter
from pandas import ExcelFile
from keras.utils import plot_model
from keras.callbacks import CSVLogger, ModelCheckpoint
from keras.callbacks import LearningRateScheduler
import math
from keras.utils import to_categorical
from keras.layers import Dropout

from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix,accuracy_score
from numpy import savetxt



def st_dec(epoch):
    # in_lr=0.8
    # drop=0.03
    # epoch_drop=3
    # lrate=in_lr*math.pow(drop,math.floor((1+epoch)/epoch_drop))
    # if lrate< 0.001:
    #     lrate=0.001

    if epoch>50 and epoch <70:
        lr=0.01
    else:
        lr=0.0001
    return lr


def makecategorical(y,ismatlab):
    mat = np.zeros((y.shape[0], y.shape[1] * 3))
    for i in range(0, y.shape[0]):
        m = 0
        for j in range(0, y.shape[1]):
            if y[i, j] == (0+ismatlab):
                mat[i, m:(m + 3)] = [1, 0, 0]
            elif y[i, j] == (1+ismatlab):
                mat[i, m:(m + 3)] = [0, 1, 0]
            elif y[i, j] == (2+ismatlab):
                mat[i, m:(m + 3)] = [0, 0, 1]
            m = m + 3
    return mat

def fromcategorical(pred,ismatlab):
    mat = [0] * 10
    j = 0
    for i in range(0, pred.shape[1], 3):
        arg = np.argmax(pred[0, i:i + 3])
        mat[j] = arg + ismatlab
        j = j + 1
    return mat

def giveprediction(elemindex,outputdim):
    predtable=np.zeros((outputdim,1))
    realtable=np.zeros((outputdim,1))
    testvalues = np.array([x_test[elemindex]])
    for nnind in range(0, outputdim):
        predtable[nnind]=model.predict_classes(testvalues)
        realtable[nnind]= y_test[elemindex, nnind]
    # print('predictions are:', predtable.T)
    # print('real data are:  ', realtable.T)
    # print('  ')
    res = np.array_equal(predtable.T, realtable.T)
    return res


csv_logger = CSVLogger('training.log', separator=',', append=False)
##
# df = pd.read_excel('datatable class=1 N=8 atan.xlsx', sheet_name='Sheet1')
df = pd.read_excel('datatable-class0-N8-testclassify-big-15240-randomized.xlsx')
ismatlab=1
conv1d_filter=0
nndata=np.array(df)
# for sf in range(0,360):
#     np.random.shuffle(nndata)


testdata_num=765
inputdim=8
offset=0 # only if inputdim is different from 8. then offset is (8-inputdim)
outputdim=8
x=nndata[0:(nndata.shape[0]-testdata_num),0:inputdim]
x_test=nndata[(nndata.shape[0]-testdata_num):nndata.shape[0],0:inputdim]

y          =np.zeros( (nndata.shape[0]-testdata_num,outputdim) )
y_test     =np.zeros( (testdata_num,outputdim) )
ycat       =np.zeros( (nndata.shape[0]-testdata_num,3,outputdim) ) # !!!!! 3 for k=[-1,1], 4 for k=[-1,0,1] !!!!
ycat_test  =np.zeros( (testdata_num,3,outputdim) )  # !!!!!! 3 for k=[-1,1], 4 for k=[-1,0,1] !!!!!
for i in range(0,outputdim):
    y[:,i]=nndata[0:(nndata.shape[0]-testdata_num),i+inputdim+offset]
    y_test[:,i]=nndata[(nndata.shape[0]-testdata_num):nndata.shape[0],i+inputdim+offset]
    # ycat[:,:,i]=to_categorical(y[:,i])
    # ycat_test[:, :, i] = to_categorical(y_test[:, i])
    ycat[:,:,i]=label_binarize(y[:,i],classes=[1,2,3])
    ycat_test[:, :, i] = label_binarize(y_test[:, i],classes=[1,2,3])

# ycat=to_categorical(y)
# ycat_test=to_categorical(y_test)


scaler = StandardScaler()
scaler.fit(x)
x=scaler.transform(x)
scaler.fit(x_test)
x_test=scaler.transform(x_test)

if conv1d_filter==1:
    xnew=np.reshape(x,(x.shape[0],x.shape[1],1))
    xnew_test=np.reshape(x_test,(x_test.shape[0],x_test.shape[1],1))
else:
    xnew = x
    xnew_test = x_test

model=[0]*outputdim
history=[0]*outputdim
mult=1

plt.figure(figsize=(2.5*1,3.5*8))
plt.legend(['train loss', 'val loss'], loc=9, borderaxespad=-2.5, ncol=2)
for nnind in range(0,outputdim):
    ker_in='glorot_uniform'
    model[nnind] = Sequential()

    if conv1d_filter==1:
        model[nnind].add(Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(xnew.shape[1], xnew.shape[2]), padding='same'))
        model[nnind].add(BatchNormalization())
        model[nnind].add(Conv1D(filters=64, kernel_size=3, activation='relu', padding='same'))
        model[nnind].add(MaxPooling1D(pool_size=2))
        model[nnind].add(Flatten())




    model[nnind].add(Dense(64,activation='relu',input_dim=x.shape[1])) #the first line of code that adds the first Dense layer is doing 2 things, defining the input or visible layer (dim=8) and the first hidden layer (dim=12)
    # model[nnind].add(LeakyReLU(alpha=.001))
    model[nnind].add(BatchNormalization())
    # model[nnind].add(Dropout(0.8))
    model[nnind].add(Dense(64,activation='relu'))
    # model[nnind].add(LeakyReLU(alpha=.001))
    model[nnind].add(BatchNormalization())
    model[nnind].add(Dense(64,activation='relu'))
    # model[nnind].add(LeakyReLU(alpha=.001))
    model[nnind].add(BatchNormalization())


    hh=model[nnind].add(Dense(ycat[:,:,nnind].shape[1], activation='sigmoid'))

    # sgd = optimizers.SGD(lr=1e-5,momentum=0.1,clipnorm=1,nesterov=True)
    rmsprop=optimizers.RMSprop(lr=0.0001, rho=0.9)
    model[nnind].compile(loss='categorical_crossentropy', optimizer=rmsprop, metrics=['mae','accuracy'])
    lr = LearningRateScheduler(st_dec,1)
    # mcp_save = ModelCheckpoint('bestmodel_'+str(nnind)+'.hdf5', save_best_only=True, monitor='val_acc', mode='max')
    history[nnind]=model[nnind].fit(xnew, ycat[:,:,nnind], epochs=100, validation_data=[xnew_test,ycat_test[:,:,nnind]], batch_size=64, shuffle=True, callbacks=[csv_logger],verbose=1)

    plt.subplot(outputdim,1,nnind+1)
    plt.plot(history[nnind].history['loss'],'--k')
    plt.plot(history[nnind].history['val_loss'],'k')
    plt.show()
    plt.pause(0.0001)
    print(np.min(history[nnind].history['loss']))
    score = model[nnind].evaluate(xnew_test, ycat_test[:, :, nnind], verbose=0)
    print('Test accuracy for', nnind, 'NN :', score[2])
    mult = mult * score[2]

    #--------- SAVE MODEL -----------
    model_json = model[nnind].to_json()
    with open("model "+str(nnind)+" .json", "w") as json_file:
        json_file.write(model_json)
    model[nnind].save_weights("model "+str(nnind)+" .h5")
    print("Saved model to disk")
    #--------------------------------


for nnind in range(0, outputdim):
    savetxt('/Users/charalamposp/Documents/Δικά μου - Κώδικες/KERAS/10NNLosses/'+str(nnind)+'model_loss.csv', history[nnind].history['loss'], delimiter=',')
    savetxt('/Users/charalamposp/Documents/Δικά μου - Κώδικες/KERAS/10NNLosses/'+str(nnind)+'model_valloss.csv', history[nnind].history['val_loss'], delimiter=',')


##
#---------------- LOAD MODEL -------------------
mult=1
loaded_model=[0]*outputdim
for nnind in range(0,outputdim):
    json_file = open('model '+str(nnind)+' .json', 'r')
    loaded_model_json = json_file.read()
    loaded_model[nnind] = model_from_json(loaded_model_json)
    loaded_model[nnind].load_weights("model "+str(nnind)+" .h5")

    loaded_model[nnind].compile(loss='categorical_crossentropy', optimizer='rmsprop', metrics=['mae','accuracy'])
    score = loaded_model[nnind].evaluate(xnew_test, ycat_test[:, :, nnind], verbose=0)
    print('Test accuracy for '+str(nnind)+' NN :', score[2])
    mult=mult*score[2]
json_file.close()
print('the product of accuracies is: '+str(mult))
##
print('the product of accuracies is: '+str(100*mult)+'%')
for nnind in range(0, outputdim):
    score = model[nnind].evaluate(xnew_test, ycat_test[:, :, nnind], verbose=0)
    print('Test accuracy for', nnind, 'NN :', score[2])

##-------------------------------------------------
# EVALUATE TOTAL
truecounter=0
predd=np.array([0]*testdata_num)
predictionall=np.zeros(y_test.shape)


for i in range (0,testdata_num):
    testvalues = np.array([xnew_test[i]])
    predictiontot=np.array([0]*outputdim)
    predd[i] =model[nnind].predict_classes(testvalues)
    for nnind in range(0, outputdim):
        predictiontot[nnind]=model[nnind].predict_classes(testvalues)

    realoutput=y_test[i]
    predictionall[i,:]=predictiontot+1
    res=np.allclose(predictiontot+1, realoutput.T, rtol=1e-02, atol=1e-02)
    if res==1:
        truecounter=truecounter+1
        print(realoutput.T)

print('Total accuracy:',100*truecounter/testdata_num,'%')

CM = np.zeros((3,3,8))
for i in range(0,8):
    CM[:,:,i]=confusion_matrix(y_test[:,i],predictionall[:,i])
    print(accuracy_score(y_test[:,i],predictionall[:,i]))
savetxt('/Users/charalamposp/Documents/Δικά μου - Κώδικες/KERAS/CM8NN.csv', np.reshape(CM,3*3*8,1), delimiter=',')

# predict the null motion
selectedtask=np.array([[0 , 0, 0, 0, 0.6178, -0.7863, 0, 0]])
predictednullmotion=np.zeros((outputdim,1))
for nnind in range(0, outputdim):
    predictednullmotion[nnind] = model[nnind].predict_classes(selectedtask)+1

##
# for nnind in range(0,outputdim):
#     loadedmodel = models.load_model('bestmodel_'+str(nnind)+'.hdf5')
#     score = loadedmodel.evaluate(xnew_test, ycat_test[:, :, 0], verbose=0)
#     print('Test accuracy for', nnind, 'NN :', score[2])



# elemindex=0
# truecounter=0
# for elemindex in range(0,100):
#     res=giveprediction(elemindex,outputdim)
#     if res==True:
#         truecounter=truecounter+1
# print(truecounter/100)




#evaluate the model
# for nnind in range(0,outputdim):


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################



############### SAVE MODEL #############
# serialize model to JSON
model_json = model.to_json()
with open("model.json", "w") as json_file:
    json_file.write(model_json)
# serialize weights to HDF5
model.save_weights("model.h5")
print("Saved model to disk")
########################################


######## LOAD MODEL ########
log_data = pd.read_csv('training.log', sep=',', engine='python')
plt.plot(log_data.loss)
plt.plot(log_data.mean_absolute_error)
# load json and create model
json_file = open('model.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
# load weights into new model
loaded_model.load_weights("model.h5")


print("Loaded model from disk")
# evaluate loaded model on test data
loaded_model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mae'])

# make class predictions with the model
testvalues=np.array([x[1]])
# testvalues=np.array([[0.7071, -0.7071,0,0]])
predictions = loaded_model.predict(testvalues)
print('predictions are:')
predictions2=[0]*10
for i in range(0,predictions.shape[1]):
    if predictions[0,i]<=np.min(predictions)+(np.max(predictions)-np.min(predictions))/3:
        predictions2[i]=0
    elif predictions[0,i]<=np.min(predictions)+2*(np.max(predictions)-np.min(predictions))/3:
        predictions2[i] = 1
    else:
        predictions2[i] = 2

print(np.array(predictions2))

print('real data are:')
print(y[1])


##### ALL #####
testvalues=x
# testvalues=np.array([[0.7071, -0.7071,0,0]])
pred = loaded_model.predict(testvalues)
modelpred=np.zeros((pred.shape[0], predictions.shape[0]))
for j in range(0,pred.shape[0]):
    predictions2=[0]*10
    predictions=pred[j,:]
    for i in range(0,predictions.shape[0]):
        if predictions[i]<=np.min(predictions)+(np.max(predictions)-np.min(predictions))/3:
            predictions2[i]=0
        elif predictions[i]<=np.min(predictions)+2*(np.max(predictions)-np.min(predictions))/3:
            predictions2[i] = 1
        else:
            predictions2[i] = 2

    modelpred[j,:]=np.array(predictions2)

print('predictions are:')
print(modelpred)
print('real data are:')
print(y)



