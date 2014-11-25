3 : 0''
if. UNAME-:'Linux' do.
  LIBCLUST=:  'libcluster.so'
elseif. UNAME-:'Darwin' do.
  'platform not supported' 13!:8[10
elseif. UNAME-: 'Win' do.
  'platform not supported' 13!:8[10
elseif. do.
  'platform not supported' 13!:8[10
end.
)
