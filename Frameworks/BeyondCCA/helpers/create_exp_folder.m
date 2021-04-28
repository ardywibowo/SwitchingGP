function exp_name = create_exp_folder(exp_name, exp_root)
  exp_nmb = 1;
  while exist(strcat(exp_root,exp_name,'_',num2str(exp_nmb)), 'dir') == 7
    exp_nmb = exp_nmb + 1;
  end
  exp_name = strcat( exp_name, '_', num2str(exp_nmb) );
  mkdir( strcat( exp_root, exp_name ) )
end