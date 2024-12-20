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
    target_size=(256, 256),
    batch_size=12,
    class_mode='categorical'
)

validation_generator = validation_datagen.flow_from_directory(
    'C://Users//Lenovo//Desktop//archive//CUB_200_2011//dataset//test',
    target_size=(256, 256),
    batch_size=12,
    class_mode='categorical'
)

# 加载预训练模型
final_model_path = os.path.join(model_save_dir, 'final_model_savedmodel')
model2 = tf.keras.models.load_model(final_model_path)

# 定义学习率调度器
initial_lr = 0.0001  # 减小初始学习率以微调模型
decay_rate = 0.9
decay_step = 2

def lr_scheduler(epoch, lr):
    if epoch % decay_step == 0 and epoch > 0:
        return lr * decay_rate
    return lr

# 定义优化器和学习率衰减回调函数
optimizer = tf.keras.optimizers.legacy.SGD(learning_rate=initial_lr, momentum=0.9, decay=1e-6)
lr_decay_callback = tf.keras.callbacks.LearningRateScheduler(lr_scheduler)

# 编译模型
model2.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])

# 继续训练模型
history = model2.fit(
    train_generator,
    validation_data=validation_generator,
    epochs=30,  # 可根据需要调整训练轮数
    callbacks=[lr_decay_callback]
)

# 最终保存整个模型为 SavedModel 格式
final_model_path = os.path.join(model_save_dir, 'final_model_savedmodel_finetuned')
model2.save(final_model_path)



