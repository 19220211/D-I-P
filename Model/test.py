import os
import tensorflow as tf
from tensorflow.keras.preprocessing import image
from PIL import Image
import matplotlib.pyplot as plt

def load_model(model_path):
    model = tf.keras.models.load_model(model_path)
    return model

def classify_images(model, image_folder, img_size=(224, 224)):
    transform = tf.keras.Sequential([
        tf.keras.layers.Resizing(img_size[0], img_size[1]),
        tf.keras.layers.Lambda(lambda x: x / 255.0)  # 归一化到 [0, 1]
    ])
    
    image_paths = [os.path.join(image_folder, f) for f in os.listdir(image_folder) if f.endswith(('.png', '.jpg', '.jpeg'))]
    results = []
    
    for img_path in image_paths:
        img = Image.open(img_path).convert('RGB')
        img_array = tf.keras.preprocessing.image.img_to_array(img)
        img_tensor = transform(img_array)[tf.newaxis, ...]  # 添加 batch 维度
        
        logits = model.predict(img_tensor)
        PRED = tf.argmax(logits, axis=1).numpy()[0]
        print("PRED:", PRED)
        
        class_map = load_class_mapping(classes_file_path)
        pred = PRED + 1  # 预测编号和 classes.txt 相对应
        preds = class_map[pred]
        print("preds:", preds)
        results.append((img_path, preds))
    
    print("results:", results)
    return results

def visualize_results(results):
    for img_path, pred in results:
        img = Image.open(img_path)
        plt.imshow(img)
        plt.title(f'Predicted Class: {pred}')
        plt.axis('off')
        plt.show()  # 可视化分类情况，不存储

def load_class_mapping(file_path):
    class_mapping = {}
    with open(file_path, 'r') as file:
        for line in file:
            index, name = line.strip().split(' ', 1)
            class_mapping[int(index)] = name
    return class_mapping

if __name__ == '__main__':
    classes_file_path = 'C://Users//Lenovo//Desktop//archive//CUB_200_2011//classes.txt'
    model_path = 'D://Model//final_model_savedmodel_finetuned'  # 替换为你的模型路径
    image_folder = 'C://Users//Lenovo//Desktop//archive//CUB_200_2011//dataset//ForTest'  # 替换为待测试的文件夹路径
    
    model = load_model(model_path)
    results = classify_images(model, image_folder)
    visualize_results(results)



