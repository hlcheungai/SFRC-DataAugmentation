% Define a new neural network for transfer learning and select layer to
% freeze
%
% Use a same name as in the pre-trained layer if the weight and bias should 
% be copied from pre-trained network.
% pre-trained layers = [ ...
%     sequenceInputLayer(13,'Name','input','Normalization','zscore')
%     gruLayer(500,'Name','gru_1','OutputMode','sequence')
%     gruLayer(500,'Name','gru_2','OutputMode','sequence')
%     gruLayer(500,'Name','gru_3','OutputMode','sequence')
%     dropoutLayer(0.5,'Name','dropout')
%     fullyConnectedLayer(6,'Name','fc')
%     regressionLayer('Name','output')
%     ];

% Define neural network
layers = [ ...
    sequenceInputLayer(13,'Name','input','Normalization','zscore')
    gruLayer(500,'Name','gru_1_new','OutputMode','sequence')
    gruLayer(500,'Name','gru_2_new','OutputMode','sequence')
    gruLayer(500,'Name','gru_3_new','OutputMode','sequence')
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(6,'Name','fc_new')
    regressionLayer('Name','output')
    ];

% Define the layers that should be freezed during training
freezeLayers = {};