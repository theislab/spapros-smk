import time
from pathlib import Path
import argparse
import pandas as pd
import scanpy as sc
import spapros as sp


# Get arguments
def get_args():
    """Get input arguments
    """
    parser = argparse.ArgumentParser(description="Run evaluation")
    
    parser.add_argument('-e', '--eval_step', help='Evaluation step (shared, pre, main, summary)', required=True, type=str)
    parser.add_argument('-m', '--metric', help='Evaluation metric', required=False, type=str)
    parser.add_argument('-d', '--data', help='Input h5ad file name', required=True, type=str)
    parser.add_argument('-r', '--results_dir', help='Output directory', required=True, type=str)
    parser.add_argument('-s', '--selection', help='Selection name', required=False, type=str)
    parser.add_argument('-E', '--evaluation_overview_file', help='Evaluation overview csv path', required=False, type=str)
    parser.add_argument('-D', '--dataset_params_file', help='Dataset parameter overview csv path', required=True, type=str)
    
    return parser.parse_args()
    
    
def main():
    """
    """
    
    # Get input arguments
    args = get_args()
    
    # Set variables
    eval_step = args.eval_step
    data_path = args.data
    eval_dataset = data_path.split("/")[-1].split(".")[0]
    output_dir = args.results_dir
    df_ds = pd.read_csv(args.dataset_params_file, index_col=0)
    ct_key = df_ds.loc[int(eval_dataset.split("_")[-1]), "ct_key"]
    if eval_step == "summary":
        df_eval = pd.read_csv(args.evaluation_overview_file, index_col=0)    
    else:
        metric = args.metric
        selection_name = args.selection    
    
    
    print("################")
    print("## EVALUATION ##")
    print("################")
    print(f"evaluating on {data_path}")
    
    
    # Load dataset
    adata = sc.read(data_path)
    print(adata)
    
    
    # Init evaluator
    start = time.time()
    print(f"Starting Evaluator setup at {start}, i.e. {time.ctime()}")
    
    evaluator = sp.ev.ProbesetEvaluator(
        adata, verbosity=0, 
        results_dir=Path(output_dir, "evaluation"), celltype_key=ct_key,
        scheme="custom", 
        metrics=[metric] if (eval_step != "summary") else [],
        reference_name=eval_dataset 
    )
    
    end = time.time()
    duration = end - start    
    print(f"Finishing Evaluator setup at {end}, i.e. {time.ctime()}")
    print(f"Evaluator setup took {duration} sec.")
    
    
    # Run evaluation
    print("\n Evaluation step: ", eval_step, "\n#######################################")
    start = time.time()
    print(f"Start evaluation at {start}, i.e. {time.ctime()}")
    
    if eval_step == "shared":
        evaluator.compute_or_load_shared_results()
    
    elif eval_step in ["pre", "main"]:
        df = pd.read_csv(Path(output_dir, "selection", f"{selection_name}.csv"), index_col=0)
        gene_set = df[df.iloc[:,0]].index.to_list()
        evaluator.evaluate_probeset(gene_set, selection_name, pre_only=(eval_step == "pre"), update_summary=False)

    elif eval_step == "summary":
        # Problem: evaluator.summary_statistics assumed that all probesets are evaluated for the same set of metrics, 
        #          which is not necessarily given here.
        # Solution: group probesets per set of metrics and set evaluator.metrics for each group separately.
        
        # Group probeset_ids based on set of metrics
        df_eval = df_eval.loc[
            (df_eval["eval_dataset"].astype(str) + "_" + df_eval["eval_data_id"].astype(str)) == eval_dataset
        ]
        probeset_ids = df_eval["selection_name"].unique().tolist()
        #probeset_ids = df_eval.loc[
        #    (df_eval["eval_dataset"].astype(str) + "_" + df_eval["eval_data_id"].astype(str)) == eval_dataset, "selection_name"
        #].unique().tolist()
        
        probeset_id_to_metrics = {
            p_id:frozenset(df_eval.loc[df_eval["selection_name"] == p_id,"metric"].values.tolist()) for p_id in probeset_ids
        }
        sets_of_metric_sets = set([metrics for metrics in probeset_id_to_metrics.values()])
        metric_sets_to_ids = {
            metrics:[p for p,ms in probeset_id_to_metrics.items() if ms == metrics] for metrics in sets_of_metric_sets
        }
        
        for metrics, set_ids in metric_sets_to_ids.items():
            # Trick to let the evaluator only load results for the current set of metrics
            evaluator.metrics = list(metrics) 
            evaluator.results = {m:{} for m in metrics}

            evaluator.summary_statistics(set_ids=set_ids)

    end = time.time()
    duration = end - start
    print(f"Finishing evaluation at {end}, i.e. {time.ctime()}")
    print(f"Evaluation took {duration} sec.")
    

if __name__ == "__main__":

    main()