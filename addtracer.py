from samples import samples

def addtracer(tracer, tlabel=''):
  for k, v in samples[tracer].items():
    gg = globals()

    gg.update({k+'{}'.format(tlabel): v})

if __name__ == '__main__':
  tracer = 'GRUSH24'
  
  addtracer(tracer, 1)

  gg     = globals()
  
  print(gg)
