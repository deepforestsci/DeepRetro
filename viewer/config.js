const config = {
  instances: [
    // Vm 2
    {
      url: 'http://localhost:5000',
      defaults: {
        model_type: 'claude4',
        advanced_prompt: true,
        model_version: 'USPTO',
        stability_flag: true,
        hallucination_check: true,
        hallucination_method: 'rule_based'
      }
    },
    // // Vm 4
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'claude37',
    //     advanced_prompt: true,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: true
    //   }
    // },
    // // VM 5
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'claude37',
    //     advanced_prompt: true,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: true
    //   }
    // },
    // // VM 3
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'deepseek',
    //     advanced_prompt: false,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: false
    //   }
    // },
    // // {
    // //   url: 'dummy ec2 url',
    // //   defaults: {
    // //     model_type: 'deepseek',
    // //     advanced_prompt: false,
    // //     model_version: 'Pistachio_100+',
    // //     stability_flag: true,
    // //     hallucination_check: true
    // //   }
    // // },
    // // VM 6
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'deepseek',
    //     advanced_prompt: false,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: true
    //   }
    // },
    // // Vm 7
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'claude3',
    //     advanced_prompt: true,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: true
    //   }
    // },
    // // Vm 10
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'claude3',
    //     advanced_prompt: true,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: true
    //   }
    // },
    // // Vm 8
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'deepseek',
    //     advanced_prompt: false,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: false
    //   }
    // },
    // // VM 9
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'deepseek',
    //     advanced_prompt: false,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: true,
    //     hallucination_check: false
    //   }
    // },
    // // VM 11
    // {
    //   url: 'dummy ec2 url',
    //   defaults: {
    //     model_type: 'claude37',
    //     advanced_prompt: true,
    //     model_version: 'Pistachio_100+',
    //     stability_flag: false,
    //     hallucination_check: true
    //   }
    // },
  ],
  endpoints: {
    retrosynthesis: '/api/retrosynthesis',
    rerun: '/api/rerun_retrosynthesis',
    partial_rerun: '/api/partial_rerun',
    save_edited: '/api/save_edited_result'
  }
};