function Indexer = imakeIndexer(Condition,animals,Cond)

Indexer = struct;

for iC = 1:length(Condition) %input from script (length of condition: 'Pre' 'Combo' 'Post1'=3) %i2 is condition we are in
    for iA = 1:length(animals) %input from groups (in workspace) %i3 is animal we are on
        counter = 0;
        counter = counter + length(Cond.(Condition{iC}){iA}); %Cond has all conditions recorded, condition(current cond type selection for current animal selection)
        % counter = the number of entries in current condition for the current animal
        if iA == 1 %for the first animal
            Indexer(1).(Condition{iC}) = counter; %the indexer structure is the size of the current counter
        elseif Indexer.(Condition{iC}) < counter %if the next animal has a structure larger than the current counter, it will make the indexer that size
            Indexer(1).(Condition{iC}) = counter;
        end
        %An indexer is now made for the i2th condition that is the size of the largest amongst the animals
    end
    % Indexer now contains all conditions at their greatest sizes
end

count = 1; %new variable 'count'
for iC = 1:length(Condition) %i2 is still which condition we are in and it's reset to 1
    if iC == 1 %first condition type
        Indexer(2).(Condition{iC}) =1; %adding a second field to the condition i2
    else
        Indexer(2).(Condition{iC}) = count;
    end
    count = count + Indexer(1).(Condition{iC}); %current count (e.g. 1)+size of condition i2 in index(e.g. 3) = (e.g. 4)
end
