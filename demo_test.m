% DEMO_TEST.m

addpath(genpath('/home/vision/wren/caffe-dilate'));
caffe.reset_all();
clc; clear;

hazy_path = ('./inputs/');
type = '*.png';
hazy_data = dir(fullfile(hazy_path,type));

model_path = './models/';

solver_file = fullfile(model_path, 'dehaze_solver_test.prototxt');
save_file = fullfile(model_path, 'dehaze_ps128_bs1.mat');

Solver = modelconfig_test( solver_file, save_file);
% Solver.iter
Solver.batchsize = 1;
for ii = 45:length(hazy_data)
    disp(ii)
    clear batch;
    hazyimg = im2double(imread(fullfile(hazy_path,hazy_data(ii).name)));

    [row, col, cha] = size(hazyimg);

% resize the input to a multiple of 8
    height = ceil(row/8)*8;
    width = ceil(col/8)*8;
    hazyimg = imresize(hazyimg, [height, width]);


    hazy_wb = RealGWbal(uint8(255*hazyimg));
    hazy_wb = hazy_wb/255;
    hazy_cont = (2*(0.5+mean(hazyimg(:)))).*(hazyimg-mean(hazyimg(:)));
    hazy_gamma = hazyimg.^2.5;

    batch(:,:,1:3) = hazyimg;
    batch(:,:,4:6) = hazy_wb;
    batch(:,:,7:9) = hazy_cont;
    batch(:,:,10:12) = hazy_gamma;

%    figure(1);subplot(221), imshow(hazyimg);subplot(222),imshow(hazy_wb);
%              subplot(223), imshow(hazy_cont); subplot(224),imshow(hazy_gamma);

    batchc = {single(batch)};

    Solver.Solver_.net.blobs('data').reshape([height, width, 12, 1]);

    tic
    activec = Solver.Solver_.net.forward(batchc);
    toc
    dehazed= activec{3};

% 	dehazed = imresize(dehazed, [size(TEST_IMAGE, 1), size(TEST_IMAGE, 2)]);
	
	imwrite(dehazed, strcat('./results/',hazy_data(ii).name(1:end-4),'_dehazed.png'));
end
%imshow(out);title('Original/Proposed/GT');
%weights=Solver.Solver_.net.get_weights();
%data=Solver.Solver_.net.get_data();
%diff=Solver.Solver_.net.get_diff();
