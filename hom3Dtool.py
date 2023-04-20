def KUBC(x, i, j, ud, bbox):
    values = np.zeros(x.shape);
    if i==j:
        values[i,:] = ud*(x[i] - bbox[i])/(bbox[i+3] - bbox[i]);
    else:
        values[i,:] = 0.5*ud*(x[j] - bbox[j])/(bbox[j+3] - bbox[j]);
        values[j,:] = 0.5*ud*(x[i] - bbox[i])/(bbox[i+3] - bbox[i]);
        
        
    return values;