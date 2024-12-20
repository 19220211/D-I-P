from PIL import Image
import tensorflow as tf
from keras.preprocessing import image
from keras.callbacks import LearningRateScheduler
import matplotlib.pyplot as plt
import os

# 设置模型保存路径并创建目录（如果不存在）
model_save_dir = 'D://Model'
os.makedirs(model_save_dir, exist_ok=True)

# 检查 GPU 是否被识别，并限制 GPU 内存增长
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(e)

# 数据增强与加载训练/验证数据集
train_datagen = tf.keras.preprocessing.image.ImageDataGenerator(
    rescale=1./255,
    rotation_range=20,
    width_shift_range=0.1,
    height_shift_range=0.1,
    shear_range=0.2,
    zoom_range=0.2,
    horizontal_flip=True,
    fill_mode='nearest'
)

validation_datagen = tf.keras.preprocessing.image.ImageDataGenerator(rescale=1./255)

train_generator = train_datagen.flow_from_directory(
    'C://Users//Lenovo//Desktop//archive//CUB_200_2011//dataset//train',
    target_size=(224, 224),
    batch_size=12,
    class_mode='categorical'
)

validation_generator = validation_datagen.flow_from_directory(
    'C://Users//Lenovo//Desktop//archive//CUB_200_2011//dataset//test',
    target_size=(224, 224),
    batch_size=12,
    class_mode='categorical'
)

# 显示训练数据集中的几张图像
train_images, train_labels = next(train_generator)
plt.figure(figsize=(10, 10))
for i in range(4):
    plt.subplot(2, 2, i+1)
    plt.imshow(train_images[i])
    plt.axis('off')
plt.show()

# 构建模型
base_model = tf.keras.applications.InceptionV3(weights='imagenet', include_top=False, input_shape=(224, 224, 3))

# 添加自定义层
x = base_model.output
x = tf.keras.layers.GlobalAveragePooling2D()(x)
x = tf.keras.layers.Dense(2048, activation='relu')(x)
x = tf.keras.layers.Dropout(0.3)(x)  # 添加Dropout层防止过拟合
x = tf.keras.layers.Dense(512, activation='relu')(x)
predictions = tf.keras.layers.Dense(200, activation='softmax')(x)
model2 = tf.keras.Model(inputs=base_model.input, outputs=predictions)

# 定义学习率调度器
initial_lr = 0.001
decay_rate = 0.9
decay_step = 2

def lr_scheduler(epoch, lr):
    if epoch % decay_step == 0 and epoch > 0:
        return lr * decay_rate
    return lr

# 定义优化器和学习率衰减回调函数
optimizer = tf.keras.optimizers.legacy.SGD(learning_rate=initial_lr, momentum=0.9, decay=1e-5)
lr_decay_callback = tf.keras.callbacks.LearningRateScheduler(lr_scheduler)

# 编译模型
model2.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])

# 训练模型
history = model2.fit(
    train_generator,
    validation_data=validation_generator,
    epochs=70,
    callbacks=[lr_decay_callback]
)

# 最终保存整个模型为 SavedModel 格式
final_model_path = os.path.join(model_save_dir, 'final_model_savedmodel')
model2.save(final_model_path)