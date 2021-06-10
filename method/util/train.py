

# transform
transform = transform_dict[opt.transform]

# data loader
data_loader = DataLoader(
    Data2D(opt, transform=transform), batch_size=opt.batch_size, 
    shuffle=True, drop_last=True, num_workers=num_workers)
    
if opt.mode == 'meta':
    model = create_model(opt.model, n_cls, opt.dataset) # repeat
    ckpt = torch.load(opt.path_s)
    model.load_state_dict(ckpt['model'])

    if torch.cuda.is_available(): # repeat
        model = model.cuda()
        cudnn.benchmark = True

# model
if opt.mode != 'meta':
    if opt.mode == 'distill':
        model_t = load_teacher(opt.path_t, n_cls, opt.dataset)
        data = torch.randn(2, 3, 84, 84)
        model_t.eval()
        model.eval()
        feat_t, _ = model_t(data, is_feat=True)
        feat, _ = model(data, is_feat=True)

        module_list = nn.ModuleList([])
        module_list.append(model)
        trainable_list = nn.ModuleList([])
        trainable_list.append(model)

        criterion_cls = nn.CrossEntropyLoss()
        criterion_div = DistillKL(opt.kd_T)
        if opt.distill == 'kd':
            criterion_kd = DistillKL(opt.kd_T)
        elif opt.distill == 'contrast':
            criterion_kd = NCELoss(opt, n_data)
            embed = Embed(feat[-1].shape[1], opt.feat_dim)
            embed_t = Embed(feat_t[-1].shape[1], opt.feat_dim)
            module_list.append(embed)
            module_list.append(embed_t)
            trainable_list.append(embed)
            trainable_list.append(embed_t)
        elif opt.distill == 'attention':
            criterion_kd = Attention()
        elif opt.distill == 'hint':
            criterion_kd = HintLoss()
        else:
            raise NotImplementedError(opt.distill)

        criterion_list = nn.ModuleList([])
        criterion_list.append(criterion_cls)    # classification loss
        criterion_list.append(criterion_div)    # KL divergence loss, original knowledge distillation
        criterion_list.append(criterion_kd)     # other knowledge distillation loss
    else:
        criterion = nn.CrossEntropyLoss()

    # optimizer
    mpar = trainable_list.parameters() if opt.mode == 'distill' else model.parameters()
    if opt.adam:
        optimizer = torch.optim.Adam(mpar,
                                    lr=opt.learning_rate,
                                    weight_decay=0.0005)
    else:
        optimizer = optim.SGD(mpar,
                            lr=opt.learning_rate,
                            momentum=opt.momentum,
                            weight_decay=opt.weight_decay)

    # append teacher after optimizer to avoid weight_decay
    if opt.mode == 'distill':
        module_list.append(model_t)

    if torch.cuda.is_available():
        if opt.mode == 'distill':
            module_list.cuda()
            criterion_list.cuda()
        else:
            if opt.n_gpu > 1:
                model = nn.DataParallel(model)
            model = model.cuda()
            criterion.cuda()
        
        cudnn.benchmark = True

    if opt.mode == 'distill':
        # validate teacher accuracy
        teacher_acc, _, _ = validate(val_loader, model_t, criterion_cls, opt)
        print('teacher accuracy: ', teacher_acc)

    # tensorboard
    logger = tb_logger.Logger(logdir=opt.tb_folder, flush_secs=2)

    # set cosine annealing scheduler
    if opt.cosine:
        eta_min = opt.learning_rate * (opt.lr_decay_rate ** 3)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, opt.epochs, eta_min, -1)

# train
if opt.mode != 'meta':

    epoch_ = 0
    save_file = os.path.join(opt.save_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch_).zfill(3)))
    save_checkpoint(model, optimizer, save_file, epoch_)

    ckpts = os.listdir(opt.save_folder)
    ckpts.sort()
    print(ckpts)
    if 'ckpt' in ckpts[-1]:
        model, optimizer, epoch_ = load_checkpoint(model, optimizer, os.path.join(opt.save_folder, ckpts[-1]))
    print(epoch_)

    for epoch in range(epoch_ + 1, opt.epochs + 1):
        
        if opt.cosine:
            scheduler.step()
        else:
            adjust_learning_rate(epoch, opt, optimizer)
        print("==> training...")

        time1 = time.time()
        if opt.mode == 'distill':
            train_acc, train_loss = train(epoch, train_loader, module_list, criterion_list, optimizer, opt)
        else:
            train_acc, train_loss = train(epoch, train_loader, model, criterion, optimizer, opt)
        
        time2 = time.time()
        print('epoch {}, total time {:.2f}'.format(epoch, time2 - time1))

        logger.log_value('train_acc', train_acc, epoch)
        logger.log_value('train_loss', train_loss, epoch)

        if opt.mode == 'distill':
            test_acc, test_acc_top5, test_loss = validate(val_loader, model, criterion_cls, opt)
        else:
            test_acc, test_acc_top5, test_loss = validate(val_loader, model, criterion, opt)

        logger.log_value('test_acc', test_acc, epoch)
        logger.log_value('test_acc_top5', test_acc_top5, epoch)
        logger.log_value('test_loss', test_loss, epoch)

        # regular saving
        if epoch % opt.save_freq == 0:
            print('==> Saving...')
            save_file = os.path.join(opt.save_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch).zfill(3)))
            save_checkpoint(model, optimizer, save_file, epoch, opt.n_gpu)

    # save the last model
    save_file = os.path.join(opt.save_folder, '{}_last.pth'.format(opt.model))
    save_checkpoint(model, optimizer, save_file, opt.epochs, opt.n_gpu)

# meta
if opt.mode == 'meta':
    # evalation
    start = time.time()
    val_acc, val_std = meta_test(model, meta_valloader)
    val_time = time.time() - start
    print('val_acc: {:.4f}, val_std: {:.4f}, time: {:.1f}'.format(val_acc, val_std,
                                                                  val_time))

    start = time.time()
    val_acc_feat, val_std_feat = meta_test(model, meta_valloader, use_logit=False)
    val_time = time.time() - start
    print('val_acc_feat: {:.4f}, val_std: {:.4f}, time: {:.1f}'.format(val_acc_feat,
                                                                       val_std_feat,
                                                                       val_time))

    start = time.time()
    test_acc, test_std = meta_test(model, meta_testloader)
    test_time = time.time() - start
    print('test_acc: {:.4f}, test_std: {:.4f}, time: {:.1f}'.format(test_acc, test_std,
                                                                    test_time))

    start = time.time()
    test_acc_feat, test_std_feat = meta_test(model, meta_testloader, use_logit=False)
    test_time = time.time() - start
    print('test_acc_feat: {:.4f}, test_std: {:.4f}, time: {:.1f}'.format(test_acc_feat,
                                                                         test_std_feat,
                                                                         test_time))