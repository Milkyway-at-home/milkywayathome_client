#import "MilkywayViewController.h"

#import <AudioToolbox/AudioToolbox.h>

#import <sys/types.h>
#import <sys/sysctl.h>
#import <stdio.h>

@implementation MilkywayViewController

@synthesize progressView, outputView, statusLabel, submitButton, startButton;

- (void)dealloc
{
    [super dealloc];
}

#pragma mark - View lifecycle

- (void)viewDidLoad
{
    outputView.text = @"";

    NSString * cachesDirectory = [NSSearchPathForDirectoriesInDomains(NSCachesDirectory, NSUserDomainMask, YES) objectAtIndex:0];
    chdir([cachesDirectory UTF8String]);
    logPath = [[cachesDirectory stringByAppendingPathComponent:@"milkyway.log"] retain];

    freopen([logPath UTF8String], "w+", stderr);
}

- (void)viewDidUnload
{
    [super viewDidUnload];
}

- (BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)interfaceOrientation
{
    return YES;
}

#pragma mark - Machine ID

- (NSString *)platform
{
    size_t size;
    sysctlbyname("hw.machine", NULL, &size, NULL, 0);
    char * machine = malloc(size);
    sysctlbyname("hw.machine", machine, &size, NULL, 0);
    NSString * platform = [NSString stringWithCString:machine encoding:NSUTF8StringEncoding];
    free(machine);
    return platform;
}

#pragma mark - Callbacks

- (void)failedToFetch
{
    statusLabel.text = @"Failed to fetch workunit...";
    submitButton.hidden = YES;
    startButton.hidden = NO;

    [progressTimer invalidate];
    progressTimer = nil;
}

- (void)runMilkywayAtHome:(id)obj
{
    NSAutoreleasePool * pool = [[NSAutoreleasePool alloc] init];

    // THIS IS TERRIFYING

    NSString * wuPath = @"http://milkyway.cs.rpi.edu/milkyway/download/milkyway_iphone/iphone_workunits/";

    NSString * jobIndex = [NSString stringWithContentsOfURL:[NSURL URLWithString:[wuPath stringByAppendingString:@"units.txt"]]
                                                   encoding:NSUTF8StringEncoding
                                                      error:nil];

    if(!jobIndex || [jobIndex isEqualToString:@""])
    {
        [self failedToFetch];

        return;
    }

    NSString * cachesDirectory = [NSSearchPathForDirectoriesInDomains(NSCachesDirectory, NSUserDomainMask, YES) objectAtIndex:0];

    NSArray * jobs = [jobIndex componentsSeparatedByCharactersInSet:[NSCharacterSet newlineCharacterSet]];
    NSString * job = [jobs objectAtIndex:random()%[jobs count]];

    NSArray * jobParts = [job componentsSeparatedByString:@" | "];
    NSArray * firstArgs = [[jobParts objectAtIndex:0] componentsSeparatedByCharactersInSet:[NSCharacterSet whitespaceCharacterSet]];
    NSArray * secondArgs = [[jobParts objectAtIndex:1] componentsSeparatedByCharactersInSet:[NSCharacterSet whitespaceCharacterSet]];

    NSURL * firstFileRemoteName = [NSURL URLWithString:[wuPath stringByAppendingString:[firstArgs objectAtIndex:0]]];
    NSURL * secondFileRemoteName = [NSURL URLWithString:[wuPath stringByAppendingString:[firstArgs objectAtIndex:1]]];
    NSString * npValue = [firstArgs objectAtIndex:2];

    NSString * firstFileContents = [NSString stringWithContentsOfURL:firstFileRemoteName encoding:NSUTF8StringEncoding error:nil];
    NSString * secondFileContents = [NSString stringWithContentsOfURL:secondFileRemoteName encoding:NSUTF8StringEncoding error:nil];

    if(!firstFileContents || [firstFileContents isEqualToString:@""] || !secondFileContents || [secondFileContents isEqualToString:@""])
    {
        [self failedToFetch];

        return;
    }

    [firstFileContents writeToFile:[cachesDirectory stringByAppendingPathComponent:@"mw-firstfile.txt"] atomically:YES encoding:NSUTF8StringEncoding error:nil];
    [secondFileContents writeToFile:[cachesDirectory stringByAppendingPathComponent:@"mw-secondfile.txt"] atomically:YES encoding:NSUTF8StringEncoding error:nil];

    NSMutableArray * arguments = [[NSMutableArray alloc] init];
    [arguments addObject:@"milkyway"];
    [arguments addObject:@"-c"];
    [arguments addObject:@"-a"];
    [arguments addObject:[cachesDirectory stringByAppendingPathComponent:@"mw-firstfile.txt"]];
    [arguments addObject:@"-s"];
    [arguments addObject:[cachesDirectory stringByAppendingPathComponent:@"mw-secondfile.txt"]];
    [arguments addObject:@"-np"];
    [arguments addObject:npValue];
    [arguments addObject:@"-p"];
    [arguments addObjectsFromArray:secondArgs];

    const char ** argv = (const char **)calloc([arguments count], sizeof(const char *));

    int i = 0;

    for(NSString * arg in arguments)
    {
        argv[i] = [arg UTF8String];
        i++;
    }

    _iphone_main([arguments count], argv);

    [pool drain];
}

- (void)threadFinished
{
    statusLabel.text = @"Finished! Submit the result.";

    submitButton.hidden = NO;
    startButton.hidden = NO;
    [progressTimer invalidate];
    progressTimer = nil;

    AudioServicesPlaySystemSound(kSystemSoundID_Vibrate);
}

- (void)updateProgress:(NSTimer *)timer
{
    if(timer != progressTimer)
        return;

    progressView.progress = _milkywaySeparationGlobalProgress;
    statusLabel.text = [NSString stringWithFormat:@"Running (%.4f%%)...", progressView.progress * 100.0, nil];

    fflush(stderr);

    outputView.text = [NSString stringWithFormat:@"Platform: %@\n%@",
                       [self platform],
                       [NSString stringWithContentsOfFile:logPath encoding:NSUTF8StringEncoding error:nil]];

    if(progressView.progress >= 1.0f)
    {
        statusLabel.text = @"Almost done!";
    }

    if(![milkywayThread isExecuting])
    {
        [self threadFinished];
    }
}

- (IBAction)start:(id)sender
{
    milkywayThread = [[NSThread alloc] initWithTarget:self selector:@selector(runMilkywayAtHome:) object:nil];
    [milkywayThread start];

    statusLabel.text = @"Running...";
    outputView.text = @"";
    submitButton.hidden = YES;
    startButton.hidden = YES;
    progressView.progress = 0.0f;

    truncate([logPath UTF8String], 0);

    progressTimer = [NSTimer scheduledTimerWithTimeInterval:1.0f target:self selector:@selector(updateProgress:) userInfo:nil repeats:YES];
}

- (IBAction)submit:(id)sender
{
    if([MFMailComposeViewController canSendMail])
    {
        MFMailComposeViewController * mailViewController = [[MFMailComposeViewController alloc] init];
        mailViewController.mailComposeDelegate = self;
        [mailViewController setToRecipients:[NSArray arrayWithObject:@"milkyway4ios@gmail.com"]];
        [mailViewController setSubject:@"milkyway@home iOS Result"];
        [mailViewController setMessageBody:[outputView text] isHTML:NO];

        [self presentModalViewController:mailViewController animated:YES];
        [mailViewController release];
    }
}

-(void)mailComposeController:(MFMailComposeViewController*)controller didFinishWithResult:(MFMailComposeResult)result error:(NSError *)error
{
    [self dismissModalViewControllerAnimated:YES];
}

@end
