import torch
import json
import esm
import pandas as pd

import argparse
from fairscale.nn.data_parallel import FullyShardedDataParallel as FSDP
from fairscale.nn.wrap import enable_wrap, wrap

def main(args):
    data_ = json.load(open(args.json_file, 'r'))
    
    sequences = []
    data_lst = []
    n = 0
    for data in data_:
        protein = data['Sequence']
        n += 1
        if len(protein) < 1000:
            sequence = (str(n), protein)
#           print(sequence)
            sequences.append(sequence)
            data_lst.append(data)

    print(len(sequences))
    
    with open('./KM_combination_1000_length_esm2_t33_650M.json', 'w') as outfile :
        json.dump(data_lst, outfile, indent=4)
    
    # load esm2 model
    # initialize the model with FSDP wrapper
    # init the distributed world with world_size 1
    url = "tcp://localhost:12349"
    torch.distributed.init_process_group(backend="nccl", init_method=url, world_size=1, rank=0)

    # download model data from the hub
    model_name = "esm2_t33_650M_UR50D"     
    model_data, regression_data = esm.pretrained._download_model_and_regression_data(model_name)

    # initialize the model with FSDP wrapper
    fsdp_params = dict(
        mixed_precision=True,
        flatten_parameters=True,
        state_dict_device=torch.device("cpu"),  # reduce GPU mem usage
        cpu_offload=True,  # enable cpu offloading
    )
    with enable_wrap(wrapper_cls=FSDP, **fsdp_params):
        model, vocab = esm.pretrained.load_model_and_alphabet_core(
            model_name, model_data, regression_data
        )
        batch_converter = vocab.get_batch_converter()
        model.eval()

        # Wrap each layer in FSDP separately
        for name, child in model.named_children():
            if name == "layers":
                for layer_name, layer in child.named_children():
                    wrapped_layer = wrap(layer)
                    setattr(child, layer_name, wrapped_layer)
        model = wrap(model)
    
    sequence_representations = []
    m = 0
    for seq in sequences:
        m += 1
        print(m, len(seq[1]))
        
        batch_labels, batch_strs, batch_tokens = batch_converter([seq])
        batch_tokens = batch_tokens.cuda()
        batch_lens = (batch_tokens != vocab.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33]
        
        for i, tokens_len in enumerate(batch_lens):
            sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))
    
    seq_embedding = pd.DataFrame(sequence_representations, dtype = float)
    seq_embedding.to_csv(args.output, sep='\t')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--json_file', type=str, help = 'Required! json file containing protein sequences.')
    parser.add_argument('-o', '--output', type=str, help = 'Required! output file')
    args = parser.parse_args()
    main(args)
