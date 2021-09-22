import os
import argparse

from transform import transform_names
from model import model_names

def parse_options():

    parser = argparse.ArgumentParser('argument for training')

    # train/test stage
    parser.add_argument('--mode', type=str, default='pretrain', choices=['pretrain', 'distill', 'meta'])

    # model
    parser.add_argument('--model', type=str, default='setr', choices=model_names)
    parser.add_argument('--model_dir', type=str, default='./model', help='model directory')
    parser.add_argument('--model_path', type=str,  help='where to save/load model')
    parser.add_argument('--preload_model', action='store_false', 
                        help="preload model to continue training or as final model if file ends in \'_final.pth\')")
    parser.add_argument('--n_class', type=int, default=5, help='max number of classes')
    parser.add_argument('--depth', type=int, default=5, help='encoder depth')
    parser.add_argument('--dim', type=int, default=200, help='side length of square image data')

    # data
    
    parser.add_argument('--data_dir', type=str, default='./data', help='data directory')
    parser.add_argument('--x_2D', type=str, default='x_2Ddenscat, x_2Dcontour', help='delimited list of input folder names in data_dir')
    parser.add_argument('--y_2D', type=str, default='y_2D', help='output folder in data_dir')

    parser.add_argument('--transform', type=str, default='A', choices=transform_names)

    parser.add_argument('--tb_dir', type=str, default='./tensorboard', help='tensorboard directory')
    
    # optimization
    parser.add_argument('--save_freq', type=int, default=10, help='pretrain: save model every save_freq epochs')
    parser.add_argument('--num_workers', type=int, default=4, help='number of workers to use')
    parser.add_argument('--batch_size', type=int, default=64, help='pretrain: batch size')
    parser.add_argument('--epochs', type=int, default=100, help='number of training epochs')

    parser.add_argument('--learning_rate', type=float, default=0.05, help='learning rate')
    parser.add_argument('--lr_decay_epochs', type=str, default='60, 80', help='delimited list of where to decay lr')
    parser.add_argument('--lr_decay_rate', type=float, default=0.1, help='decay rate for learning rate')
    parser.add_argument('--weight_decay', type=float, default=5e-4, help='weight decay')
    parser.add_argument('--momentum', type=float, default=0.9, help='momentum')
    
    # cosine annealing
    parser.add_argument('--cosine', action='store_true', help='using cosine annealing')
    
    # meta train/test
    parser.add_argument('--data_scat', type=str, default='pregnancy/07_FoxP3CD25_CD4Tcell', help='meta: dataset/scatterplot folders')
    parser.add_argument('--shot_dir', type=str, default='./data/x_2Ddensity_euclidean_rankkmed', help='meta: directory with shot names as filenames')
    parser.add_argument('--n_shots', type=int, default=1, metavar='N',
                        help='meta: number of support samples')
    parser.add_argument('--meta_batch_size', type=int, default=1, metavar='meta_batch_size',
                        help='meta: batch size')
    
    parser.add_argument('-t', '--trial', type=str, default='1', help='experiment id')

    # distillation
    parser.add_argument('--distill', type=str, default='kd', choices=['kd', 'contrast', 'hint', 'attention'])

    parser.add_argument('-r', '--gamma', type=float, default=1, help='weight for classification')
    parser.add_argument('-a', '--alpha', type=float, default=0, help='weight balance for KD')
    parser.add_argument('-b', '--beta', type=float, default=0, help='weight balance for other losses')

    # KL distillation
    parser.add_argument('--kd_T', type=float, default=4, help='temperature for KD distillation')

    # NCE distillation
    parser.add_argument('--feat_dim', default=128, type=int, help='feature dimension')
    parser.add_argument('--nce_k', default=16384, type=int, help='number of negative samples for NCE')
    parser.add_argument('--nce_t', default=0.07, type=float, help='temperature parameter for softmax')
    parser.add_argument('--nce_m', default=0.5, type=float, help='momentum for non-parametric updates')

    # ======================================================
    opt = parser.parse_args('')

    opt.gpu = torch.cuda.is_available()
    opt.n_gpu = torch.cuda.device_count()

    opt.lr_decay_epochs = [int(dpi) for dpi in delim_pars(opt.lr_decay_epochs)]
    opt.x_2D = delim_pars(opt.x_2D)
    
    # dirs
    opt.model_name = '{}_depth:{}_dim:{}_epoch:{}_trans:{}_lr:{}_decay:{}'.format(
    opt.model, opt.depth, opt.dim, opt.epochs, opt.learning_rate, opt.weight_decay, opt.transform)
    # if opt.mode == 'distill':
    #     opt.model_t = get_teacher_name(opt.path_t)
    #     opt.model_name = 's:{}_t:{}_trans:{}_d:{}_r:{}_a:{}_b:{}'.format(
    #         opt.model, opt.model_t, opt.distill, 
    #         opt.gamma, opt.alpha, opt.beta, opt.transform)
    if opt.cosine:
        opt.model_name = '{}_cosine'.format(opt.model_name)
    if opt.adam:
        opt.model_name = '{}_useAdam'.format(opt.model_name)
    opt.model_name = '{}_{}'.format(opt.model_name, opt.trial)

    if 'meta' not in opt.mode:
        opt.tb_folder = os.path.join(opt.tb_path, opt.model_name)
        os.makedirs(opt.tb_folder, exist_ok=True) # exist_ok only on python 3.2+
    else:
        opt.transform = 'B' # only aug size
        opt.batch_size = opt.meta_batch_size

        opt.model_name_meta = '{}_META_shots:{}'.format(opt.model_name, opt.n_shots)

        # shot_dir = ./data/x_2Ddensity_euclidean_rankkmed / pregnancy/07_FoxP3CD25_CD4Tcell / 1-5, 10, 15, 20
        opt.shot_dir = os.path.join(opt.shot_dir, opt.data_scat, str(opt.n_shots))
        
    opt.save_dir = opt.model_dir
    os.makedirs(opt.model_dir, exist_ok=True)

    return opt