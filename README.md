# A collection of python classes for anavar (ref)

## About


## API

### Control file classes

All control file classes require the use of ```.set_data()```

As a bare minimum the site frequency data, sample size and number of callable sites must be given.

The input site frequency data takes the form of a dictrionary where:

```python
site_frequencies = [0.9, 0.6, 0.5, 0.4, 0.4, 0.8, 0.3]
sample_size = 17
site_freq_format = {'key': (site_frequencies, sample_size)}
```

Where the number of entries and the keys in the dictionary depend on the model:

* SNP_1 keys: 'SNP'
* INDEL_1 keys: 'INS', 'DEL'
* gBGC_GLEMIN_EXTENDED_M1* keys: 'neutral_SNPs', 'ws_SNPs', 'sw_SNPs'
* neutralINDEL_vs_selectedINDEL: 'neutral_INS', 'neutral_DEL', 'selected_INS', 'selected_DEL'