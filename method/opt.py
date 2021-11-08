import os
import argparse
import torch

from models import model_names

def parse_options():
    parser = argparse.ArgumentParser()
    
    # train/test stage
    parser.add_argument('--mode', type=str, default='pretrain', choices=['pretrain', 'meta'])

    # root dir
    parser.add_argument('--root_dir', type=str, default='/home/aya43/flowMagic_data', help='root directory')
    # parser.add_argument('--root_dir', type=str, default='/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data', help='root directory')

    # model
    parser.add_argument('--model', type=str, default='deeplab3', choices=model_names)
    parser.add_argument('--model_dir', type=str, default='model', help='model directory')
    parser.add_argument('--new_weights', action='store_false', 
                        help="don't preload model to continue training or as final model if file ends in \'_final.pth\')")
    parser.add_argument('--n_class', type=int, default=6, help='max number of classes')
    parser.add_argument('--dim', type=int, default=256, help='side length of square image data')

    # data
    parser.add_argument('--data_dir', type=str, default='data/2D', help='data directory')
    parser.add_argument('--x_2D', type=str, default='x_2Ddenscat,x_2Dcontour', help=', delimited list of input (channel) folder names in data_dir')
    parser.add_argument('--y_2D', type=str, default='y_2D,x_2Ddiscrete,y_vector_', help=', delimited list of output folder names data_dir')
    parser.add_argument('--preload_data', action='store_false', help='preload all data into memory')

    parser.add_argument('--tb_dir', type=str, default='tensorboard', help='tensorboard directory')
    
    # optimization
    parser.add_argument('--save_freq', type=int, default=10, help='pretrain: save model every save_freq epochs')
    parser.add_argument('--print_freq', type=int, default=10, help='pretrain: print model score every print_freq epochs')
    
    parser.add_argument('--num_workers', type=int, default=32, help='number of workers to use')
    parser.add_argument('--batch_size', type=int, default=32, help='pretrain: batch size')
    parser.add_argument('--epochs', type=int, default=100, help='number of training epochs')
    
    parser.add_argument('--adam', action='store_true', help='using Adam GD')
    parser.add_argument('--learning_rate', type=float, default=0.001, help='learning rate')
    parser.add_argument('--lr_decay_epochs', type=str, default='60,80', help='delimited list of where to decay lr')
    parser.add_argument('--lr_decay_rate', type=float, default=0.1, help='decay rate for learning rate')
    parser.add_argument('--weight_decay', type=float, default=5e-2, help='weight decay')
    parser.add_argument('--momentum', type=float, default=0.9, help='momentum')
    
    # meta train/test
    parser.add_argument('--data_scat', type=str, default='pregnancy/07_FoxP3CD25_CD4Tcell', help='meta: dataset/scatterplot folder')
    parser.add_argument('--shot_dir', type=str, default='data/2D/x_2Ddenscat_euclidean_rankkmed', help='meta: directory with shot names as filenames')
    parser.add_argument('--n_shots', type=int, default=1, metavar='N',
                        help='meta: number of meta-training support samples')

    # ======================================================
    opt = parser.parse_args('') # opt = parser.parse_args() # from console
    opt = update_opt(opt)

    return opt

def update_opt(opt):
    opt.gpu = torch.cuda.is_available()
    opt.n_gpu = torch.cuda.device_count()
    
    if ',' in opt.lr_decay_epochs:
        opt.lr_decay_epochs = [int(dpi) for dpi in opt.lr_decay_epochs.split(',')]
    if ',' in opt.x_2D:
        opt.x_2D = opt.x_2D.split(',')
    if ',' in opt.y_2D:
        opt.y_2D = opt.y_2D.split(',')
    
    # experiment name
    opt.model_name = '{}_dim:{}_epoch:{}_lr:{}_decay:{}'.format(
    opt.model, opt.dim, opt.epochs, opt.learning_rate, opt.weight_decay)
    
    opt.adam = True  ## USE ADAM
    if opt.adam:
        opt.model_name = '{}_useAdam'.format(opt.model_name)
    
    opt.data_folder = os.path.join(opt.root_dir, opt.data_dir)
    opt.model_folder = os.path.join(opt.root_dir, opt.model_dir, opt.model_name)
    opt.tb_folder = os.path.join(opt.root_dir, opt.tb_dir, opt.model_name)
    opt.shot_folder = os.path.join(opt.root_dir, opt.shot_dir, opt.data_scat, str(opt.n_shots))
        
    os.makedirs(opt.model_folder, exist_ok=True) # exist_ok only on python 3.2+
    
    return opt