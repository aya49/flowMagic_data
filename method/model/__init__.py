from .vit_mla import model_setr

model_dict = {
    'setr': model_setr,
    'unet': None
}

model_names = list()
for name, dict_ in model_dict.items():
    model_names.append(name)