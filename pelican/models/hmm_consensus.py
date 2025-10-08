"""
Script para gerar e avaliar um HMM para filtragem de genes não consenso.
"""
import numpy as np
from scipy.ndimage import uniform_filter1d
from hmmlearn import hmm

# Função para construir trilhas de sinais
def build_signal_tracks(region_sets, genome_length, smooth_window=5):
    starts = np.zeros(genome_length)
    overlaps = np.zeros(genome_length)
    ends = np.zeros(genome_length)
    for regions in region_sets:
        for start, end in regions:
            # Ajusta índices para evitar IndexError
            start_idx = max(0, start)
            end_idx = min(genome_length - 1, end - 1)
            starts[start_idx] += 1
            ends[end_idx] += 1
            overlaps[start_idx:end_idx+1] += 1
    starts = uniform_filter1d(starts, size=smooth_window)
    ends = uniform_filter1d(ends, size=smooth_window)
    return starts, overlaps, ends

# Função para normalizar scores (softmax)
def softmax(scores):
    exp_scores = np.exp(scores - np.max(scores))
    return exp_scores / exp_scores.sum()

# Função para avaliar genes não consenso com HMM

def calculate_transition_matrix(states_seq, n_states=3):
    """
    Calcula a matriz de transição baseada na sequência de estados estimada.
    """
    transmat = np.zeros((n_states, n_states))
    for (i, j) in zip(states_seq[:-1], states_seq[1:]):
        transmat[i, j] += 1
    # Normaliza para probabilidades
    transmat = np.where(transmat.sum(axis=1, keepdims=True) == 0, 1, transmat)
    transmat = transmat / transmat.sum(axis=1, keepdims=True)
    return transmat


def evaluate_genes_hmm(starts, overlaps, ends, region_sets, n_states=3, dynamic_transmat=False, prob_threshold=0.99, transmat=None):
    """
    Executa HMM e retorna genes filtrados por probabilidade, mantendo genes únicos em regiões.
    Agora com apenas três estados: início, core e fim.
    """
    X = np.column_stack([starts, overlaps, ends])
    model = hmm.GaussianHMM(n_components=n_states, covariance_type="diag")
    model.startprob_ = np.array([1.0, 0.0, 0.0])
    model.means_ = np.array([
        [np.max(starts), np.mean(overlaps), np.min(ends)],   # início
        [np.mean(starts), np.max(overlaps), np.mean(ends)],  # core
        [np.min(starts), np.mean(overlaps), np.max(ends)],   # fim
    ])
    model.covars_ = np.tile(np.var(X, axis=0), (n_states, 1))
    if transmat is not None:
        model.transmat_ = transmat
    else:
        model.transmat_ = np.array([
            [0.95, 0.05, 0.0],
            [0.0, 0.95, 0.05],
            [0.05, 0.0, 0.95],
        ])
    states = model.predict(X)
    if dynamic_transmat:
        model.transmat_ = calculate_transition_matrix(states, n_states)
        states = model.predict(X)
    score = model.score(X)
    # Obter probabilidades de cada estado para cada posição
    logprob, posteriors = model.score_samples(X)
    # Filtrar genes por probabilidade dos estados início (0) ou fim (2) > 0.5
    filtered_genes = []
    window = 13
    prob_threshold_state = 0.6  # threshold para início/fim ESTAVA EM 0.6
    # Para garantir que genes em regiões únicas não tenham qualquer sobreposição
    all_gene_coords = []
    for regions in region_sets:
        for start, end in regions:
            all_gene_coords.append((start, end))

    def is_unique_region(start, end, all_coords):
        for s2, e2 in all_coords:
            if (s2, e2) == (start, end):
                continue
            # Se houver qualquer sobreposição
            if max(start, s2) < min(end, e2):
                return False
        return True

    for i, regions in enumerate(region_sets):
        for start, end in regions:
            region_probs = posteriors[start:end]  # shape: (gene_len, n_states)
            gene_len = end - start
            # Garante que o gene não tem qualquer sobreposição
            if is_unique_region(start, end, all_gene_coords):
                filtered_genes.append((start, end, True))
            else:
                # Probabilidades nas janelas
                initial_probs = region_probs[:window, 0]  # início
                final_probs = region_probs[-window:, 2]   # fim
                high_start = np.any(initial_probs > prob_threshold_state)
                high_end = np.any(final_probs > prob_threshold_state)
                # Gene é aceito se início ou fim tem probabilidade > 0.5
                if high_start or high_end:
                    filtered_genes.append((start, end, False))
    return filtered_genes, score, model.transmat_

if __name__ == "__main__":
    # Exemplo de uso
    # region_sets = [[(10, 50), (100, 150)], [(20, 60), (120, 170)]]
    # genome_length = 200
    # starts, overlaps, ends = build_signal_tracks(region_sets, genome_length)
    # states, score = evaluate_genes_hmm(starts, overlaps, ends)
    # print("Estados:", states)
    # print("Score:", score)
    pass
