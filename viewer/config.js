const config = {
    instances: [
        'http://ec2-18-220-15-234.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-18-118-82-88.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-18-222-202-30.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-18-118-34-4.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-3-20-205-97.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-3-135-238-102.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-18-224-139-148.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-18-117-169-71.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-3-138-103-124.us-east-2.compute.amazonaws.com:5000',
        'http://ec2-13-59-117-221.us-east-2.compute.amazonaws.com:5000'
    ],
    endpoints: {
        retrosynthesis: '/api/retrosynthesis',
        rerun: '/api/rerun_retrosynthesis',
        partial_rerun: '/api/partial_rerun'
    }
};